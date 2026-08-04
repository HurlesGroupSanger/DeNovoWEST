#!/usr/bin/env python
"""
Generate a self-contained HTML report from DeNovoWEST pipeline results.

Reads combined_tests_results.tsv and any available per-gene stats files from
the results directory, then writes a single HTML file with a sortable gene
table and expandable per-gene detail panels.
"""

####################
# Imports          #
####################

from __future__ import annotations

import json
import logging
from datetime import date
from pathlib import Path

import click
import gffutils
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from denovowest.utils.log import init_log

####################
# Constants        #
####################

COMBINED_TESTS_RELPATH = "combined_tests/combined_tests_results.tsv"
RATES_STATS_RELPATH = "rates/stats/rates_stats.json"
FILTERED_RATES_STATS_RELPATH = "rates/stats/filtered_rates_stats.json"

PVAL_COLUMNS = [
    "min_pval",
    "combined_mis_pval",
    "enrichment_all_pval",
    "enrichment_mis_pval",
    "clustering_pval",
]

# Fixed nominal threshold; the corrected threshold is supplied by the user at runtime.
NOMINAL_THRESHOLD = 5e-2

# Caution thresholds — genes exceeding either are flagged in the report.
CAUTION_ANNOTATION_THRESHOLD = 0.10
CAUTION_FILTERED_THRESHOLD = 0.10

# Report format version — bump when the HTML layout or metric definitions change.
REPORT_VERSION = "0.1.0"

# Columns always present in simulation log observed_dnms tables; anything else is a score.
_FIXED_OBSERVED_COLS: frozenset[str] = frozenset({"chrom", "pos", "ref", "alt", "consequence"})

# G2P confidence levels in descending order (highest evidence first).
CONFIDENCE_RANK: dict[str, int] = {
    "definitive": 0,
    "strong": 1,
    "moderate": 2,
    "limited": 3,
    "disputed": 4,
    "refuted": 5,
    "absent": 6,
}

# Allelic requirements that indicate dominant (heterozygous) disease — used to
# prefer monoallelic entries when a gene appears in G2P under multiple modes.
MONOALLELIC_REQUIREMENTS: frozenset[str] = frozenset(
    {
        # Current G2P format
        "monoallelic_autosomal",
        "monoallelic_X_hemizygous",
        "monoallelic_X_heterozygous",
        "monoallelic_PAR",
        # Legacy G2P format
        "monoallelic",
    }
)

####################
# Data loading     #
####################


def load_combined_tests(results_dir: Path) -> pd.DataFrame:
    """Load the combined tests results TSV.

    Args:
        results_dir: pipeline output directory.

    Returns:
        DataFrame with one row per gene.

    Raises:
        FileNotFoundError: if the expected TSV is absent.
    """

    path = results_dir / COMBINED_TESTS_RELPATH
    if not path.exists():
        raise FileNotFoundError(f"Combined tests results not found: {path}")
    return pd.read_csv(path, sep="\t")


def _detect_score_column(columns: list[str]) -> str | None:
    """Return the score column from an observed_dnms column list.

    The score column is the only one not among the fixed set of
    chrom / pos / ref / alt / consequence.

    Args:
        columns: column list from an observed_dnms split-orient dict.

    Returns:
        Column name of the score, or None if all columns are fixed.
    """

    extras = [c for c in columns if c not in _FIXED_OBSERVED_COLS]
    return extras[0] if extras else None


def _detect_predictor_name(simulation_log: dict | None) -> str | None:
    """Return the effect predictor column name from a simulation log.

    Scans the first entry that contains an observed_dnms block and extracts
    the non-fixed column name, which is the name of the computational effect
    predictor used during the simulation (e.g. 'CADD', 'AlphaMissense').

    Args:
        simulation_log: dict mapping gene_id to simulation stats, as returned
            by load_simulation_log() or load_mis_simulation_log().

    Returns:
        Predictor column name, or None if the log is absent or has no score column.
    """

    if not simulation_log:
        return None
    for entry in simulation_log.values():
        raw = (entry or {}).get("observed_dnms")
        if raw and "columns" in raw:
            name = _detect_score_column(raw["columns"])
            if name:
                return name
    return None


def _load_simulation_log_by_runtype(results_dir: Path, runtype: str) -> dict | None:
    """Load a specific enrichment simulation log by runtype subfolder name.

    Args:
        results_dir: pipeline output directory.
        runtype: name of the subfolder under ``simulation/`` (e.g. ``all-coding``
            or ``mis``).

    Returns:
        Dict mapping versioned gene_id to simulation stats, or None if absent.
    """

    path = results_dir / "simulation" / runtype / "simulation_logs.json"
    if not path.exists():
        return None
    with path.open() as fh:
        return json.load(fh)


def load_simulation_log(results_dir: Path) -> dict | None:
    """Load the all-coding enrichment simulation log.

    Args:
        results_dir: pipeline output directory.

    Returns:
        Dict mapping versioned gene_id to simulation stats, or None if absent.
    """

    return _load_simulation_log_by_runtype(results_dir, "all-coding")


def load_mis_simulation_log(results_dir: Path) -> dict | None:
    """Load the missense enrichment simulation log.

    Args:
        results_dir: pipeline output directory.

    Returns:
        Dict mapping versioned gene_id to simulation stats, or None if absent.
    """

    return _load_simulation_log_by_runtype(results_dir, "mis")


def load_rates_stats(results_dir: Path) -> dict | None:
    """Load the per-gene per-consequence annotation statistics JSON.

    This file is produced by the rates_stats Nextflow process and contains
    counts of missing scores and mutation counts per gene, consequence, and
    predictor column.

    Args:
        results_dir: pipeline output directory.

    Returns:
        Dict mapping gene_id to per-consequence annotation stats, or None if absent.
    """

    path = results_dir / RATES_STATS_RELPATH
    if not path.exists():
        return None
    with path.open() as fh:
        return json.load(fh)


def _compute_missing_pct(
    rates_stats: dict,
    gene_id: str,
    predictor: str,
    consequences: list[str] | None = None,
) -> float | None:
    """Compute the fraction of variants missing a predictor score for a gene.

    Args:
        rates_stats: full rates_stats dict as returned by load_rates_stats().
        gene_id: versioned or unversioned Ensembl gene ID.
        predictor: predictor column name to look up (e.g. ``"score"`` or
            ``"AlphaMissense_score"``).
        consequences: if provided, restrict the calculation to these consequences;
            if None, all consequences present in the gene entry are used.

    Returns:
        Fraction missing (0–1), or None if the gene or predictor is absent.
    """

    entry = rates_stats.get(gene_id)
    if entry is None:
        # Fallback: try unversioned ID
        base_id = gene_id.split(".")[0]
        entry = next((v for k, v in rates_stats.items() if k.split(".")[0] == base_id), None)
    if entry is None:
        return None

    total_mutations = 0
    total_missing = 0
    for csq, csq_stats in entry.items():
        if consequences is not None and csq not in consequences:
            continue
        if not isinstance(csq_stats, dict):
            continue
        predictor_stats = csq_stats.get(predictor)
        if not isinstance(predictor_stats, dict):
            continue
        total_mutations += csq_stats.get("nb_mutations", 0)
        total_missing += predictor_stats.get("nb_missing_score", 0)

    if total_mutations == 0:
        return None
    return total_missing / total_mutations


def load_filtered_rates_stats(results_dir: Path) -> dict | None:
    """Load the per-gene site-count statistics before and after exclusion-region filtering.

    Args:
        results_dir: pipeline output directory.

    Returns:
        Dict with keys ``nb_loci_before_filtering`` and ``nb_loci_filtered``,
        each mapping gene_id to an integer site count, or None if absent.
    """

    path = results_dir / FILTERED_RATES_STATS_RELPATH
    if not path.exists():
        return None
    with path.open() as fh:
        return json.load(fh)


def _compute_filtered_pct(filtered_rates_stats: dict, gene_id: str) -> float | None:
    """Compute the fraction of sites removed by exclusion-region filtering for a gene.

    Args:
        filtered_rates_stats: dict from load_filtered_rates_stats().
        gene_id: versioned or unversioned Ensembl gene ID.

    Returns:
        Fraction of sites filtered (0–1), or None if the gene is absent or
        the before count is zero.
    """

    before_map: dict = filtered_rates_stats.get("nb_loci_before_filtering", {})
    filtered_map: dict = filtered_rates_stats.get("nb_loci_filtered", {})

    nb_before = before_map.get(gene_id)
    nb_filtered = filtered_map.get(gene_id)

    # Fallback to unversioned ID when not found under the versioned key
    if nb_before is None:
        base_id = gene_id.split(".")[0]
        nb_before = next((v for k, v in before_map.items() if k.split(".")[0] == base_id), None)
        nb_filtered = next((v for k, v in filtered_map.items() if k.split(".")[0] == base_id), None)

    if nb_before is None or nb_before == 0:
        return None
    return (nb_filtered or 0) / nb_before


def _best_confidence(group: pd.DataFrame) -> pd.Series:
    """Return the best confidence and allelic requirement for one gene.

    Monoallelic entries are preferred; if none exist the best entry across all
    allelic requirements is returned instead.
    """
    mono = group[group["allelic_requirement"].isin(MONOALLELIC_REQUIREMENTS)]
    candidates = mono if not mono.empty else group
    best = candidates.sort_values("_rank").iloc[0]
    return pd.Series(
        {
            "g2p_confidence": best["g2p_confidence"],
            "g2p_allelic_requirement": best["allelic_requirement"],
        }
    )


def load_g2p(g2p_path: Path) -> pd.DataFrame:
    """Load a G2P CSV and return the best confidence level per gene symbol.

    Selection rules applied in order:
    1. If the gene has any monoallelic_autosomal entries, pick the highest
       confidence among those.
    2. Otherwise pick the highest confidence among all entries.

    Args:
        g2p_path: path to a G2P CSV file.

    Returns:
        DataFrame with columns [gene_symbol, g2p_confidence], one row per gene.
    """

    df = pd.read_csv(
        g2p_path,
        usecols=["gene symbol", "allelic requirement", "confidence"],
    )
    df = df.rename(
        columns={
            "gene symbol": "gene_symbol",
            "allelic requirement": "allelic_requirement",
            "confidence": "g2p_confidence",
        }
    )

    # Drop rows with unrecognised confidence values (e.g. CSV parsing artefacts)
    df = df[df["g2p_confidence"].isin(CONFIDENCE_RANK)].copy()
    df["_rank"] = df["g2p_confidence"].map(CONFIDENCE_RANK)

    result = df.groupby("gene_symbol", group_keys=False).apply(_best_confidence, include_groups=False).reset_index()
    # result columns: gene_symbol, g2p_confidence, g2p_allelic_requirement
    return result


def merge_g2p(df: pd.DataFrame, g2p_df: pd.DataFrame) -> pd.DataFrame:
    """Add a g2p_confidence column to the results dataframe by merging on gene_symbol.

    Genes absent from G2P are labelled "absent".

    Args:
        df: combined tests results.
        g2p_df: output of load_g2p().

    Returns:
        Copy of df with a new g2p_confidence column.

    Raises:
        ValueError: if the results dataframe has no gene_symbol column.
    """

    if "gene_symbol" not in df.columns:
        raise ValueError(
            "Results file has no gene_symbol column — G2P annotation requires a pipeline "
            "version that outputs gene symbols."
        )
    merged = df.merge(g2p_df, on="gene_symbol", how="left")
    merged["g2p_confidence"] = merged["g2p_confidence"].fillna("absent")
    merged["g2p_allelic_requirement"] = merged["g2p_allelic_requirement"].fillna("absent")
    return merged


####################
# GFF loading      #
####################


def load_gff_gene_coords(gff_path: Path) -> pd.DataFrame:
    """Extract gene-level coordinates from a GFF3 file using gffutils.

    Args:
        gff_path: path to a GFF3 or GFF3.gz file.

    Returns:
        DataFrame with columns [gene_id, gene_name, chrom, start, end].
        gene_id is the versioned Ensembl ID (e.g. ENSG00000186092.7).
    """

    records = []
    for feature in gffutils.DataIterator(str(gff_path)):
        if feature.featuretype != "gene":
            continue
        records.append(
            {
                "gene_id": feature.attributes.get("gene_id", [""])[0],
                "gene_name": feature.attributes.get("gene_name", [""])[0],
                "chrom": feature.chrom,
                "start": feature.start,
                "end": feature.end,
            }
        )

    df = pd.DataFrame(records)
    # Store a version-stripped ID for fallback matching when the results file
    # uses unversioned Ensembl IDs.
    df["base_gene_id"] = df["gene_id"].str.split(".").str[0]
    return df


def build_ambiguous_symbol_gene_ids(gene_coords_df: pd.DataFrame) -> dict[str, list[str]]:
    """Return gene symbols that map to more than one gene ID in the GFF.

    Args:
        gene_coords_df: output of load_gff_gene_coords().

    Returns:
        Dict mapping gene_name to a sorted list of gene_ids, only for symbols
        with two or more distinct gene IDs.
    """

    grouped = gene_coords_df.groupby("gene_name")["gene_id"].apply(lambda ids: sorted(ids.unique()))
    return {name: list(ids) for name, ids in grouped.items() if len(ids) > 1}


####################
# Caution flags    #
####################


def _check_missing_pvals(row: pd.Series) -> list[str]:
    """Flag genes where one or more p-value tests produced no result."""
    # clustering_pval and combined_mis_pval are expected NA when fewer than 2
    # missense DNMs are observed — clustering requires ≥2 observations.
    # clustering_nb_mis_variants is NaN (not 0) when the gene was absent from the
    # clustering output entirely, meaning it had zero missense DNMs.
    if "clustering_nb_mis_variants" not in row.index:
        nb_mis: int | None = None  # column absent — cannot determine
    elif pd.isna(row["clustering_nb_mis_variants"]):
        nb_mis = 0  # absent from clustering file → zero missense DNMs
    else:
        nb_mis = int(row["clustering_nb_mis_variants"])
    skip: set[str] = set()
    if nb_mis is not None and nb_mis <= 1:
        skip.update({"clustering_pval", "combined_mis_pval"})
    if nb_mis is not None and nb_mis == 0:
        skip.add("enrichment_mis_pval")

    missing = [col for col in PVAL_COLUMNS if pd.isna(row.get(col)) and col not in skip]
    if missing:
        return [f"Missing p-value(s): {', '.join(missing)}"]
    return []


def _check_missing_annotations(
    gene_id: str,
    rates_stats_log: dict,
    score_predictor: str | None,
    mis_score_predictor: str | None,
    threshold: float = CAUTION_ANNOTATION_THRESHOLD,
) -> list[str]:
    """Flag genes where a predictor annotation is missing for more than ``threshold`` of sites."""

    reasons = []
    if score_predictor:
        pct = _compute_missing_pct(rates_stats_log, gene_id, score_predictor)
        if pct is not None and pct > threshold:
            reasons.append(f"{pct * 100:.1f}% of sites are missing {score_predictor} annotation")
    if mis_score_predictor:
        pct = _compute_missing_pct(rates_stats_log, gene_id, mis_score_predictor, consequences=["missense"])
        if pct is not None and pct > threshold:
            reasons.append(f"{pct * 100:.1f}% of missense sites are missing {mis_score_predictor} annotation")
    return reasons


def _check_filtered_sites(
    gene_id: str,
    filtered_rates_stats_log: dict,
    threshold: float = CAUTION_FILTERED_THRESHOLD,
) -> list[str]:
    """Flag genes where more than ``threshold`` of sites were removed by exclusion-region filtering."""

    pct = _compute_filtered_pct(filtered_rates_stats_log, gene_id)
    if pct is not None and pct > threshold:
        return [f"{pct * 100:.1f}% of sites were filtered out based on exclusion regions"]
    return []


def _check_symbol_ambiguity(
    gene_id: str,
    gene_coords_df: pd.DataFrame,
    ambiguous_symbol_gene_ids: dict[str, list[str]],
) -> list[str]:
    """Flag genes whose GFF symbol maps to more than one gene ID.

    Args:
        gene_id: versioned or unversioned Ensembl gene ID from the results file.
        gene_coords_df: output of load_gff_gene_coords().
        ambiguous_symbol_gene_ids: output of build_ambiguous_symbol_gene_ids().

    Returns:
        A one-element list with the warning message, or an empty list.
    """

    # Locate this gene in the GFF — try versioned ID first, then unversioned
    match = gene_coords_df[gene_coords_df["gene_id"] == gene_id]
    if match.empty:
        base_id = gene_id.split(".")[0]
        match = gene_coords_df[gene_coords_df["base_gene_id"] == base_id]
    if match.empty:
        return []

    gene_name = match.iloc[0]["gene_name"]
    if gene_name not in ambiguous_symbol_gene_ids:
        return []

    # Report the other IDs that share this symbol
    all_ids = ambiguous_symbol_gene_ids[gene_name]
    other_ids = [gid for gid in all_ids if gid != gene_id]
    return [f"Gene symbol '{gene_name}' also assigned to: {', '.join(other_ids)}, DNN results may be affected."]


def _check_g2p_overlap(
    row: pd.Series,
    corrected_threshold: float,
    gene_coords_df: pd.DataFrame,
    g2p_gene_set: frozenset[str],
) -> list[str]:
    """Flag significant non-G2P genes that overlap a G2P gene in the genome.

    Only applied when the gene is significant (min_pval < corrected_threshold)
    and absent from G2P.
    """

    if row.get("g2p_confidence") != "absent":
        return []
    min_pval = row.get("min_pval")
    if pd.isna(min_pval) or float(min_pval) >= corrected_threshold:
        return []

    gene_id = row["gene_id"]
    base_id = gene_id.split(".")[0]

    # Locate this gene's coordinates — try versioned ID first, then unversioned
    query = gene_coords_df[gene_coords_df["gene_id"] == gene_id]
    if query.empty:
        query = gene_coords_df[gene_coords_df["base_gene_id"] == base_id]
    if query.empty:
        return []

    chrom, start, end = query.iloc[0][["chrom", "start", "end"]]

    # Find all GFF genes on the same chromosome that overlap this locus
    same_chrom = gene_coords_df[gene_coords_df["chrom"] == chrom]
    overlapping = same_chrom[
        (same_chrom["end"] >= start) & (same_chrom["start"] <= end) & (same_chrom["gene_id"] != gene_id)
    ]

    g2p_hits = sorted(overlapping.loc[overlapping["gene_name"].isin(g2p_gene_set), "gene_name"])
    if g2p_hits:
        return [f"Overlaps G2P gene(s): {', '.join(g2p_hits)}"]
    return []


def compute_caution_flags(
    df: pd.DataFrame,
    corrected_threshold: float = NOMINAL_THRESHOLD,
    gene_coords_df: pd.DataFrame | None = None,
    g2p_gene_set: frozenset[str] | None = None,
    ambiguous_symbol_gene_ids: dict[str, list[str]] | None = None,
    rates_stats_log: dict | None = None,
    score_predictor: str | None = None,
    mis_score_predictor: str | None = None,
    filtered_rates_stats_log: dict | None = None,
) -> pd.DataFrame:
    """Add caution flag and reasons columns to the results dataframe.

    Args:
        df: combined tests results.
        corrected_threshold: significance threshold used for the G2P overlap check.
        gene_coords_df: GFF-derived gene coordinates for overlap detection (optional).
        g2p_gene_set: set of gene symbols present in G2P (optional).
        ambiguous_symbol_gene_ids: output of build_ambiguous_symbol_gene_ids(); genes
            whose GFF symbol resolves to multiple gene IDs are flagged (optional).
        rates_stats_log: per-gene annotation statistics from load_rates_stats() (optional).
        score_predictor: all-coding predictor column name; used for the annotation check.
        mis_score_predictor: missense predictor column name; used for the annotation check.
        filtered_rates_stats_log: site-filtering stats from load_filtered_rates_stats() (optional).

    Returns:
        Copy of df with added 'caution' (bool) and 'caution_reasons' (list[str]) columns.
    """

    do_overlap_check = gene_coords_df is not None and g2p_gene_set is not None
    do_symbol_check = gene_coords_df is not None and ambiguous_symbol_gene_ids is not None
    do_annotation_check = rates_stats_log is not None
    do_filter_check = filtered_rates_stats_log is not None

    reasons_list = []
    for _, row in df.iterrows():
        gene_id = row["gene_id"]
        reasons = _check_missing_pvals(row)
        if do_annotation_check:
            reasons += _check_missing_annotations(gene_id, rates_stats_log, score_predictor, mis_score_predictor)
        if do_filter_check:
            reasons += _check_filtered_sites(gene_id, filtered_rates_stats_log)
        if do_symbol_check:
            reasons += _check_symbol_ambiguity(gene_id, gene_coords_df, ambiguous_symbol_gene_ids)
        if do_overlap_check:
            reasons += _check_g2p_overlap(row, corrected_threshold, gene_coords_df, g2p_gene_set)
        reasons_list.append(reasons)

    df = df.copy()
    df["caution_reasons"] = reasons_list
    df["caution"] = df["caution_reasons"].map(lambda r: len(r) > 0)
    return df


####################
# Data serialisation
####################


def _safe_float(val) -> float | None:
    """Return float or None for NaN/missing values."""
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    return float(val)


def _safe_int(val) -> int | None:
    """Return int or None for NaN/missing values."""
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    return int(val)


def _safe_str(val) -> str | None:
    """Return str or None for NaN/missing values."""
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    return str(val)


def _build_gene_records(
    df: pd.DataFrame,
    simulation_log: dict | None = None,
    mis_simulation_log: dict | None = None,
    corrected_threshold: float = NOMINAL_THRESHOLD,
    rates_stats_log: dict | None = None,
    score_predictor: str | None = None,
    mis_score_predictor: str | None = None,
    filtered_rates_stats_log: dict | None = None,
) -> list[dict]:
    """Serialise the results dataframe into a JSON-friendly list of dicts.

    Args:
        df: combined tests results with caution columns.
        simulation_log: optional dict from load_simulation_log(), keyed by gene_id.
        mis_simulation_log: optional dict from load_mis_simulation_log(); the score
            column is extracted and joined onto each missense variant by position.
        corrected_threshold: significance threshold; the full per-variant DNM list is
            included only for significant genes that are not already strong/definitive
            monoallelic in G2P, keeping the JSON compact for well-established genes.
        rates_stats_log: optional dict from load_rates_stats(); used to compute the
            fraction of variants missing a predictor annotation per gene.
        score_predictor: column name of the all-coding effect predictor (e.g. ``"score"``).
        mis_score_predictor: column name of the missense effect predictor
            (e.g. ``"AlphaMissense_score"``).
        filtered_rates_stats_log: optional dict from load_filtered_rates_stats(); used to
            compute the fraction of sites removed by exclusion-region filtering per gene.

    Returns:
        List of dicts, one per gene, ready for json.dumps().
    """

    # Build base-id indexes for unversioned ID fallback lookups
    sim_base_idx: dict[str, dict] = {}
    if simulation_log:
        for gid, entry in simulation_log.items():
            sim_base_idx[gid.split(".")[0]] = entry

    mis_base_idx: dict[str, dict] = {}
    if mis_simulation_log:
        for gid, entry in mis_simulation_log.items():
            mis_base_idx[gid.split(".")[0]] = entry

    records = []
    for _, row in df.iterrows():
        gene_id = row["gene_id"]
        base_id = gene_id.split(".")[0]

        min_pval = _safe_float(row.get("min_pval"))
        g2p_confidence = _safe_str(row.get("g2p_confidence")) or "absent"
        g2p_allelic_req = _safe_str(row.get("g2p_allelic_requirement")) or "absent"
        # Genes already well-established as dominant (strong/definitive monoallelic) are
        # excluded — the variant list is most useful for novel or weakly supported candidates.
        is_established_dominant = (
            g2p_confidence in {"definitive", "strong"} and g2p_allelic_req in MONOALLELIC_REQUIREMENTS
        )
        is_significant = min_pval is not None and min_pval < corrected_threshold and not is_established_dominant

        # Pull simulation stats (versioned ID first, then unversioned fallback)
        sim_entry = (simulation_log or {}).get(gene_id) or sim_base_idx.get(base_id)

        # Count observed DNMs per consequence and sum scores, from the simulation log
        obs_per_csq: dict[str, int] | None = None
        obs_score_per_csq: dict[str, float] | None = None
        observed_dnms: list[dict] | None = None
        if sim_entry and "observed_dnms" in sim_entry:
            raw = sim_entry["observed_dnms"]
            try:
                obs_df = pd.DataFrame(data=raw["data"], columns=raw["columns"])
                score_col = _detect_score_column(raw["columns"])
                obs_per_csq = obs_df["consequence"].value_counts().to_dict()
                # Observed score sum per consequence (available for all genes, not just significant)
                if score_col:
                    obs_score_per_csq = obs_df.groupby("consequence")[score_col].sum().to_dict()
                # Full per-variant list only for significant genes, sorted by score descending
                if is_significant:
                    if score_col:
                        obs_df = obs_df.sort_values(score_col, ascending=False).rename(columns={score_col: "score"})
                    observed_dnms = obs_df.to_dict(orient="records")
            except (KeyError, ValueError):
                pass

        # Parse the missense simulation log entry for this gene
        mis_entry = (mis_simulation_log or {}).get(gene_id) or mis_base_idx.get(base_id)

        # Parse the observed_dnms block from the mis simulation log (independent of
        # whether the full per-variant list is included, so obs/exp scores are always available)
        mis_df: pd.DataFrame | None = None
        mis_score_col: str | None = None
        if mis_entry and "observed_dnms" in mis_entry:
            mis_raw = mis_entry["observed_dnms"]
            try:
                mis_df = pd.DataFrame(data=mis_raw["data"], columns=mis_raw["columns"])
                mis_score_col = _detect_score_column(mis_raw["columns"])
            except (KeyError, ValueError):
                pass

        # Observed missense score = sum of predictor scores across all missense DNMs
        obs_mis_score: float | None = None
        if mis_df is not None and mis_score_col:
            obs_mis_score = _safe_float(mis_df[mis_score_col].sum())

        # Expected missense score from the mis simulation log
        exp_mis_score: float | None = None
        if mis_entry:
            exp_mis_score = _safe_float((mis_entry.get("expected_score_per_consequence") or {}).get("missense"))

        # Attach per-variant missense scores to the observed DNM list (significant genes only)
        if observed_dnms and mis_df is not None and mis_score_col:
            mis_lookup: dict[tuple, float | None] = {
                (r["chrom"], r["pos"], r["ref"], r["alt"]): _safe_float(r[mis_score_col])
                for r in mis_df.to_dict(orient="records")
            }
            for v in observed_dnms:
                key = (v["chrom"], v["pos"], v["ref"], v["alt"])
                v["mis_score"] = mis_lookup.get(key)

        # Fraction of variants missing a predictor score, from the rates_stats log
        pct_missing_all: float | None = None
        pct_missing_mis: float | None = None
        if rates_stats_log is not None:
            if score_predictor:
                pct_missing_all = _compute_missing_pct(rates_stats_log, gene_id, score_predictor)
            if mis_score_predictor:
                pct_missing_mis = _compute_missing_pct(
                    rates_stats_log, gene_id, mis_score_predictor, consequences=["missense"]
                )

        # Fraction of sites removed by exclusion-region filtering
        pct_filtered: float | None = None
        if filtered_rates_stats_log is not None:
            pct_filtered = _compute_filtered_pct(filtered_rates_stats_log, gene_id)

        records.append(
            {
                "gene_id": gene_id,
                "base_id": base_id,
                "gene_symbol": _safe_str(row.get("gene_symbol")),
                "g2p_confidence": _safe_str(row.get("g2p_confidence")),
                "g2p_allelic_requirement": _safe_str(row.get("g2p_allelic_requirement")),
                "min_pval": _safe_float(row.get("min_pval")),
                "combined_mis_pval": _safe_float(row.get("combined_mis_pval")),
                "enrichment_all_pval": _safe_float(row.get("enrichment_all_pval")),
                "enrichment_mis_pval": _safe_float(row.get("enrichment_mis_pval")),
                "clustering_pval": _safe_float(row.get("clustering_pval")),
                "observed_all": _safe_float(row.get("observed_all")),
                "expected_all": _safe_float(row.get("expected_all")),
                "observed_mis": _safe_float(row.get("observed_mis")),
                "expected_mis": _safe_float(row.get("expected_mis")),
                "clustering_nb_mis_variants": _safe_int(row.get("clustering_nb_mis_variants")),
                "clustering_mode": _safe_str(row.get("clustering_mode")),
                "clustering_dist": _safe_float(row.get("clustering_dist")),
                "nb_observed_dnms": _safe_int(sim_entry["nb_observed_dnms"]) if sim_entry else None,
                "nb_expected_dnms": _safe_float(sim_entry["nb_expected_dnms"]) if sim_entry else None,
                "sim_expected_per_csq": sim_entry.get("expected_prob_per_consequence") if sim_entry else None,
                "sim_expected_score_per_csq": sim_entry.get("expected_score_per_consequence") if sim_entry else None,
                "obs_per_csq": obs_per_csq,
                "obs_score_per_csq": obs_score_per_csq,
                "obs_mis_score": obs_mis_score,
                "exp_mis_score": exp_mis_score,
                "pct_missing_all": pct_missing_all,
                "pct_missing_mis": pct_missing_mis,
                "pct_filtered": pct_filtered,
                "observed_dnms": observed_dnms,
                "caution": bool(row["caution"]),
                "caution_reasons": row["caution_reasons"],
            }
        )
    return records


####################
# HTML generation  #
####################

_TEMPLATE_DIR = Path(__file__).parent / "templates"


def build_html(
    df: pd.DataFrame,
    results_dir: Path,
    run_date: str,
    corrected_threshold: float,
    simulation_log: dict | None = None,
    mis_simulation_log: dict | None = None,
    rates_stats_log: dict | None = None,
    filtered_rates_stats_log: dict | None = None,
    score_predictor: str | None = None,
    mis_score_predictor: str | None = None,
) -> str:
    """Render the complete HTML report via the Jinja2 template.

    Args:
        df: combined tests results with caution columns.
        results_dir: pipeline output directory (used in the report title).
        run_date: ISO-format date string for the report header.
        corrected_threshold: multiple-testing corrected significance threshold.
        simulation_log: optional dict from load_simulation_log().
        mis_simulation_log: optional dict from load_mis_simulation_log().
        rates_stats_log: optional dict from load_rates_stats(); used to compute
            per-gene missing-annotation percentages.
        filtered_rates_stats_log: optional dict from load_filtered_rates_stats(); used to
            compute the fraction of sites removed by exclusion-region filtering.
        score_predictor: all-coding predictor name; detected from simulation_log when None.
        mis_score_predictor: missense predictor name; detected from mis_simulation_log when None.

    Returns:
        Self-contained HTML document as a string.

    Note:
        gene_json is passed as a plain string and rendered with the ``| safe``
        filter in the template, so autoescape does not double-encode it.
    """

    n_total = len(df)
    n_significant = int((df["min_pval"].fillna(1.0) < corrected_threshold).sum())
    n_caution = int(df["caution"].sum())

    # Fall back to detecting predictor names if not supplied by the caller
    if score_predictor is None:
        score_predictor = _detect_predictor_name(simulation_log)
    if mis_score_predictor is None:
        mis_score_predictor = _detect_predictor_name(mis_simulation_log)

    gene_records = _build_gene_records(
        df,
        simulation_log,
        mis_simulation_log,
        corrected_threshold,
        rates_stats_log=rates_stats_log,
        score_predictor=score_predictor,
        mis_score_predictor=mis_score_predictor,
        filtered_rates_stats_log=filtered_rates_stats_log,
    )
    # Escape </ to prevent a premature </script> tag inside the JSON block
    gene_json = json.dumps(gene_records).replace("</", "<\\/")

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(["html"]),
    )
    template = env.get_template("report.html.j2")
    return template.render(
        title=f"DeNovoWEST Report — {results_dir.name}",
        results_dir=str(results_dir),
        run_date=run_date,
        n_total=n_total,
        n_significant=n_significant,
        n_caution=n_caution,
        corrected_threshold=corrected_threshold,
        nominal_threshold=NOMINAL_THRESHOLD,
        report_version=REPORT_VERSION,
        gene_json=gene_json,
        score_predictor=score_predictor,
        mis_score_predictor=mis_score_predictor,
    )


####################
# CLI entry point  #
####################


@click.command()
@click.argument(
    "results_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--threshold",
    "-t",
    type=float,
    required=True,
    help="Multiple-testing corrected significance threshold (e.g. 2.5e-6 for Bonferroni).",
)
@click.option(
    "--g2p",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="G2P CSV file (e.g. G2P_DD_2026-06-11.csv). Confidence is merged on gene symbol.",
)
@click.option(
    "--gff",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "GFF3 file (plain or .gz). Used to suppress false-positive splice-site warnings "
        "for single-exon genes, and (when --g2p is also provided) to flag significant "
        "non-G2P genes that overlap a known G2P gene."
    ),
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Output HTML file. Defaults to <results_dir>/report.html.",
)
def report(results_dir: Path, threshold: float, g2p: Path | None, gff: Path | None, output: Path | None) -> None:
    """Generate an HTML report summarising DeNovoWEST pipeline results.

    RESULTS_DIR is the root output directory produced by the pipeline
    (e.g. /path/to/benchmarking_output/2026-06-10).

    Args:
        results_dir: pipeline output root.
        threshold: multiple-testing corrected significance threshold.
        g2p: optional G2P CSV for confidence annotation.
        gff: optional GFF3 for G2P overlap detection.
        output: destination HTML file (default: <results_dir>/report.html).
    """

    init_log()
    logger = logging.getLogger("logger")

    if output is None:
        output = results_dir / "report.html"

    run_date = date.today().isoformat()
    logger.info(f"Loading results from {results_dir}")
    logger.info(f"Corrected significance threshold: {threshold:.2e}")

    # Load combined test results
    df = load_combined_tests(results_dir)

    g2p_df = None
    if g2p is not None:
        try:
            g2p_df = load_g2p(g2p)
            df = merge_g2p(df, g2p_df)
            logger.info(f"Merged G2P annotations for {g2p_df['gene_symbol'].nunique()} genes from {g2p.name}")
        except ValueError as e:
            logger.warning(f"Skipping G2P annotation: {e}")

    gene_coords_df = None
    ambiguous_symbol_gene_ids = None
    if gff is not None:
        logger.info(f"Loading gene coordinates from {gff.name}")
        gene_coords_df = load_gff_gene_coords(gff)
        logger.info(f"Loaded coordinates for {len(gene_coords_df)} genes")
        ambiguous_symbol_gene_ids = build_ambiguous_symbol_gene_ids(gene_coords_df)
        if ambiguous_symbol_gene_ids:
            logger.info(f"Found {len(ambiguous_symbol_gene_ids)} gene symbol(s) mapping to multiple gene IDs in GFF")
        if g2p_df is None:
            logger.warning("--gff G2P overlap detection requires --g2p; only symbol-ambiguity checks will be applied")

    simulation_log = load_simulation_log(results_dir)
    if simulation_log is not None:
        logger.info(f"Loaded all-coding simulation log for {len(simulation_log)} genes")
    else:
        logger.warning("All-coding simulation log not found — observed/expected DNM counts will be unavailable")

    mis_simulation_log = load_mis_simulation_log(results_dir)
    if mis_simulation_log is not None:
        logger.info(f"Loaded missense simulation log for {len(mis_simulation_log)} genes")
    else:
        logger.info("Missense simulation log not found — missense scores will be omitted from variant table")

    rates_stats_log = load_rates_stats(results_dir)
    if rates_stats_log is not None:
        logger.info(f"Loaded annotation stats for {len(rates_stats_log)} genes")
    else:
        logger.info(
            "Annotation stats (rates_stats.json) not found — missing-annotation percentages will be unavailable"
        )

    filtered_rates_stats_log = load_filtered_rates_stats(results_dir)
    if filtered_rates_stats_log is not None:
        logger.info("Loaded exclusion-region filtering stats")
    else:
        logger.info(
            "Filtering stats (filtered_rates_stats.json) not found — site-filtering percentages will be unavailable"
        )

    # Detect predictor names early so they can inform both caution flags and the report
    score_predictor = _detect_predictor_name(simulation_log)
    mis_score_predictor = _detect_predictor_name(mis_simulation_log)

    # Build G2P gene set for overlap detection (only when both --g2p and --gff were provided)
    g2p_gene_set = frozenset(g2p_df["gene_symbol"]) if (g2p_df is not None and gene_coords_df is not None) else None

    # Compute caution flags and generate HTML
    df = compute_caution_flags(
        df,
        threshold,
        gene_coords_df,
        g2p_gene_set,
        ambiguous_symbol_gene_ids,
        rates_stats_log=rates_stats_log,
        score_predictor=score_predictor,
        mis_score_predictor=mis_score_predictor,
        filtered_rates_stats_log=filtered_rates_stats_log,
    )
    n_caution = int(df["caution"].sum())
    logger.info(f"{n_caution} gene(s) flagged with caution")

    html = build_html(
        df,
        results_dir,
        run_date,
        threshold,
        simulation_log,
        mis_simulation_log,
        rates_stats_log,
        filtered_rates_stats_log,
        score_predictor=score_predictor,
        mis_score_predictor=mis_score_predictor,
    )
    output.write_text(html, encoding="utf-8")
    logger.info(f"Report written to {output}")


if __name__ == "__main__":
    report()
