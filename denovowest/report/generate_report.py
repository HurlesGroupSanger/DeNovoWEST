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

import json
import logging
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Optional

import click
import gffutils
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from denovowest.utils.log import init_log

####################
# Constants        #
####################

COMBINED_TESTS_RELPATH = "combined_tests/combined_tests_results.tsv"
NB_RECORDS_RELPATH = "rates/stats/nb_records_per_gene_per_csq.tsv"
RATES_RELPATH = "rates/stats/rates_per_gene_per_csq.tsv"

PVAL_COLUMNS = [
    "min_pval",
    "combined_mis_pval",
    "enrichment_all_pval",
    "enrichment_mis_pval",
    "clustering_pval",
]

FUNCTIONAL_CONSEQUENCES = [
    "missense",
    "stop_gained",
    "start_lost",
    "stop_lost",
    "splice_region",
    "splice_acceptor",
    "splice_donor",
]

ALL_CONSEQUENCES = ["synonymous"] + FUNCTIONAL_CONSEQUENCES

# Splice consequences only exist in multi-exon genes; single-exon genes legitimately
# have zero rates/records for these and should not be flagged.
SPLICE_CONSEQUENCES: frozenset[str] = frozenset(
    {
        "splice_region",
        "splice_acceptor",
        "splice_donor",
    }
)

# Fixed nominal threshold; the corrected threshold is supplied by the user at runtime.
NOMINAL_THRESHOLD = 5e-2

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


def _load_wide_stats(path: Path) -> Optional[pd.DataFrame]:
    """Load a wide-format stats TSV (genes as rows, consequences as columns)."""
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "gene_id"
    return df


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


def load_nb_records(results_dir: Path) -> Optional[pd.DataFrame]:
    """Load the per-gene per-consequence record count matrix, if present."""
    return _load_wide_stats(results_dir / NB_RECORDS_RELPATH)


def load_rates(results_dir: Path) -> Optional[pd.DataFrame]:
    """Load the per-gene per-consequence mutation rate matrix, if present."""
    return _load_wide_stats(results_dir / RATES_RELPATH)


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

    result = df.groupby("gene_symbol", group_keys=False).apply(_best_confidence).reset_index()
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


def load_gff_single_cds_genes(gff_path: Path) -> frozenset[str]:
    """Return unversioned gene IDs for genes with only a single exon position.

    Single-exon genes have no introns and therefore no splice sites.  Zero rates
    or record counts for splice consequences are biologically expected for these
    genes and should not be treated as data-quality issues.

    Args:
        gff_path: path to a GFF3 or GFF3.gz file.

    Returns:
        frozenset of unversioned gene IDs (e.g. ``ENSG00000186092``) for genes
        whose annotated exon positions collapse to a single unique (start, end) pair
        across all transcripts.

    Note:
        This function iterates the GFF independently from load_gff_gene_coords,
        so two I/O passes occur when both are called.
    """

    exon_positions: dict[str, set] = defaultdict(set)
    for feature in gffutils.DataIterator(str(gff_path)):
        if feature.featuretype != "CDS":
            continue
        gene_id = feature.attributes.get("gene_id", [""])[0]
        if gene_id:
            base_id = gene_id.split(".")[0]
            exon_positions[base_id].add((feature.start, feature.end))

    return frozenset(gid for gid, positions in exon_positions.items() if len(positions) == 1)


####################
# Caution flags    #
####################


def _check_missing_pvals(row: pd.Series) -> list[str]:
    """Flag genes where one or more p-value tests produced no result."""
    # clustering_pval and combined_mis_pval are expected NA when fewer than 2
    # missense DNMs are observed — clustering requires ≥2 observations.
    # nb_mis_variants is NaN (not 0) when the gene was absent from the clustering
    # output entirely, meaning it had zero missense DNMs.
    if "nb_mis_variants" not in row.index:
        nb_mis: Optional[int] = None  # column absent — cannot determine
    elif pd.isna(row["nb_mis_variants"]):
        nb_mis = 0  # absent from clustering file → zero missense DNMs
    else:
        nb_mis = int(row["nb_mis_variants"])
    skip: set[str] = set()
    if nb_mis is not None and nb_mis <= 1:
        skip.update({"clustering_pval", "combined_mis_pval"})
    if nb_mis is not None and nb_mis == 0:
        skip.add("enrichment_mis_pval")

    missing = [col for col in PVAL_COLUMNS if pd.isna(row.get(col)) and col not in skip]
    if missing:
        return [f"Missing p-value(s): {', '.join(missing)}"]
    return []


def _check_zero_rates(
    gene_id: str,
    rates_df: Optional[pd.DataFrame],
    single_cds_gene_ids: Optional[frozenset[str]] = None,
) -> list[str]:
    """Flag genes with a zero total mutation rate for any functional consequence."""
    if rates_df is None or gene_id not in rates_df.index:
        return []
    row = rates_df.loc[gene_id]

    # Single-exon genes have no splice sites — suppress expected zeroes
    base_id = gene_id.split(".")[0]
    skip = SPLICE_CONSEQUENCES if (single_cds_gene_ids and base_id in single_cds_gene_ids) else frozenset()

    zero_csqs = [csq for csq in FUNCTIONAL_CONSEQUENCES if csq in row.index and row[csq] == 0.0 and csq not in skip]
    if zero_csqs:
        return [f"Zero mutation rate for: {', '.join(zero_csqs)}"]
    return []


def _check_zero_records(
    gene_id: str,
    nb_records_df: Optional[pd.DataFrame],
    single_cds_gene_ids: Optional[frozenset[str]] = None,
) -> list[str]:
    """Flag genes with no annotated records for any functional consequence."""
    if nb_records_df is None or gene_id not in nb_records_df.index:
        return []
    row = nb_records_df.loc[gene_id]

    # Single-exon genes have no splice sites — suppress expected zeroes
    base_id = gene_id.split(".")[0]
    skip = SPLICE_CONSEQUENCES if (single_cds_gene_ids and base_id in single_cds_gene_ids) else frozenset()

    zero_csqs = [csq for csq in FUNCTIONAL_CONSEQUENCES if csq in row.index and row[csq] == 0 and csq not in skip]
    if zero_csqs:
        return [f"No annotated records for: {', '.join(zero_csqs)}"]
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
    nb_records_df: Optional[pd.DataFrame],
    rates_df: Optional[pd.DataFrame],
    corrected_threshold: float = NOMINAL_THRESHOLD,
    gene_coords_df: Optional[pd.DataFrame] = None,
    g2p_gene_set: Optional[frozenset[str]] = None,
    single_cds_gene_ids: Optional[frozenset[str]] = None,
    ambiguous_symbol_gene_ids: Optional[dict[str, list[str]]] = None,
) -> pd.DataFrame:
    """Add caution flag and reasons columns to the results dataframe.

    Args:
        df: combined tests results.
        nb_records_df: per-gene per-consequence record counts (optional).
        rates_df: per-gene per-consequence mutation rates (optional).
        corrected_threshold: significance threshold used for the G2P overlap check.
        gene_coords_df: GFF-derived gene coordinates for overlap detection (optional).
        g2p_gene_set: set of gene symbols present in G2P (optional).
        single_cds_gene_ids: unversioned gene IDs for single-exon genes; splice-site
            warnings are suppressed for these genes (optional).
        ambiguous_symbol_gene_ids: output of build_ambiguous_symbol_gene_ids(); genes
            whose GFF symbol resolves to multiple gene IDs are flagged (optional).

    Returns:
        Copy of df with added 'caution' (bool) and 'caution_reasons' (list[str]) columns.
    """

    do_overlap_check = gene_coords_df is not None and g2p_gene_set is not None
    do_symbol_check = gene_coords_df is not None and ambiguous_symbol_gene_ids is not None

    reasons_list = []
    for _, row in df.iterrows():
        gene_id = row["gene_id"]
        reasons = (
            _check_missing_pvals(row)
            + _check_zero_rates(gene_id, rates_df, single_cds_gene_ids)
            + _check_zero_records(gene_id, nb_records_df, single_cds_gene_ids)
        )
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


def _safe_float(val) -> Optional[float]:
    """Return float or None for NaN/missing values."""
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    return float(val)


def _safe_int(val) -> Optional[int]:
    """Return int or None for NaN/missing values."""
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    return int(val)


def _safe_str(val) -> Optional[str]:
    """Return str or None for NaN/missing values."""
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    return str(val)


def _build_gene_records(
    df: pd.DataFrame,
    nb_records_df: Optional[pd.DataFrame],
    rates_df: Optional[pd.DataFrame],
) -> list[dict]:
    """Serialise the results dataframe into a JSON-friendly list of dicts.

    Args:
        df: combined tests results with caution columns.
        nb_records_df: optional per-consequence record counts.
        rates_df: optional per-consequence mutation rates.

    Returns:
        List of dicts, one per gene, ready for json.dumps().
    """

    records = []
    for _, row in df.iterrows():
        gene_id = row["gene_id"]
        base_id = gene_id.split(".")[0]

        # Build per-consequence breakdown for detail panel
        csq_data = {}
        for csq in ALL_CONSEQUENCES:
            entry = {}
            if nb_records_df is not None and gene_id in nb_records_df.index and csq in nb_records_df.columns:
                v = nb_records_df.loc[gene_id, csq]
                entry["nb_records"] = _safe_int(v)
            if rates_df is not None and gene_id in rates_df.index and csq in rates_df.columns:
                v = rates_df.loc[gene_id, csq]
                entry["rate"] = _safe_float(v)
            if entry:
                csq_data[csq] = entry

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
                "nb_mis_variants": _safe_int(row.get("nb_mis_variants")),
                "clustering_dist": _safe_float(row.get("clustering_dist")),
                "caution": bool(row["caution"]),
                "caution_reasons": row["caution_reasons"],
                "csq_data": csq_data,
            }
        )
    return records


####################
# HTML generation  #
####################

_TEMPLATE_DIR = Path(__file__).parent / "templates"


def build_html(
    df: pd.DataFrame,
    nb_records_df: Optional[pd.DataFrame],
    rates_df: Optional[pd.DataFrame],
    results_dir: Path,
    run_date: str,
    corrected_threshold: float,
) -> str:
    """Render the complete HTML report via the Jinja2 template.

    Args:
        df: combined tests results with caution columns.
        nb_records_df: optional per-consequence record counts.
        rates_df: optional per-consequence mutation rates.
        results_dir: pipeline output directory (used in the report title).
        run_date: ISO-format date string for the report header.
        corrected_threshold: multiple-testing corrected significance threshold.

    Returns:
        Self-contained HTML document as a string.

    Note:
        gene_json is passed as a plain string and rendered with the ``| safe``
        filter in the template, so autoescape does not double-encode it.
    """

    n_total = len(df)
    n_significant = int((df["min_pval"].fillna(1.0) < corrected_threshold).sum())
    n_caution = int(df["caution"].sum())

    gene_records = _build_gene_records(df, nb_records_df, rates_df)
    # Escape </ to prevent a premature </script> tag inside the JSON block
    gene_json = json.dumps(gene_records, indent=2).replace("</", "<\\/")

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
        gene_json=gene_json,
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
def report(
    results_dir: Path, threshold: float, g2p: Optional[Path], gff: Optional[Path], output: Optional[Path]
) -> None:
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

    # Load all available data sources
    df = load_combined_tests(results_dir)
    nb_records_df = load_nb_records(results_dir)
    rates_df = load_rates(results_dir)

    g2p_df = None
    if g2p is not None:
        try:
            g2p_df = load_g2p(g2p)
            df = merge_g2p(df, g2p_df)
            logger.info(f"Merged G2P annotations for {g2p_df['gene_symbol'].nunique()} genes from {g2p.name}")
        except ValueError as e:
            logger.warning(f"Skipping G2P annotation: {e}")

    gene_coords_df = None
    single_cds_gene_ids = None
    ambiguous_symbol_gene_ids = None
    if gff is not None:
        logger.info(f"Loading gene coordinates from {gff.name}")
        gene_coords_df = load_gff_gene_coords(gff)
        logger.info(f"Loaded coordinates for {len(gene_coords_df)} genes")
        single_cds_gene_ids = load_gff_single_cds_genes(gff)
        logger.info(
            f"Identified {len(single_cds_gene_ids)} single-exon gene(s) — splice-site warnings suppressed for these"
        )
        ambiguous_symbol_gene_ids = build_ambiguous_symbol_gene_ids(gene_coords_df)
        if ambiguous_symbol_gene_ids:
            logger.info(f"Found {len(ambiguous_symbol_gene_ids)} gene symbol(s) mapping to multiple gene IDs in GFF")
        if g2p_df is None:
            logger.warning(
                "--gff G2P overlap detection requires --g2p; only single-exon splice suppression will be applied"
            )

    if nb_records_df is not None:
        logger.info(f"Loaded record counts for {len(nb_records_df)} genes")
    else:
        logger.warning("Record count stats not found — caution flags will be limited")

    if rates_df is not None:
        logger.info(f"Loaded mutation rates for {len(rates_df)} genes")
    else:
        logger.warning("Mutation rate stats not found — caution flags will be limited")

    # Build G2P gene set for overlap detection (only when both --g2p and --gff were provided)
    g2p_gene_set = frozenset(g2p_df["gene_symbol"]) if (g2p_df is not None and gene_coords_df is not None) else None

    # Compute caution flags and generate HTML
    df = compute_caution_flags(
        df,
        nb_records_df,
        rates_df,
        threshold,
        gene_coords_df,
        g2p_gene_set,
        single_cds_gene_ids,
        ambiguous_symbol_gene_ids,
    )
    n_caution = int(df["caution"].sum())
    logger.info(f"{n_caution} gene(s) flagged with caution")

    html = build_html(df, nb_records_df, rates_df, results_dir, run_date, threshold)
    output.write_text(html, encoding="utf-8")
    logger.info(f"Report written to {output}")


if __name__ == "__main__":
    report()
