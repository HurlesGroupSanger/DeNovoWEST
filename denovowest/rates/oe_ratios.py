#!/usr/bin/env python

import json
import logging
from pathlib import Path
from typing import Optional

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import poisson

from denovowest.utils.helpers import compute_expected_number_of_mutations_stats, ensembl_ids_map_version
from denovowest.utils.io_helpers import load_gff
from denovowest.utils.log import init_log
from denovowest.utils.params import ConsequenceGroups

##############
# CONSTANTS  #
##############

CONSEQUENCES = ConsequenceGroups.OE_RATIO

# Label used for genes present in the rates matrix but absent from the user-supplied
# gene-category file; they are always shown as a separate category in the plots.
_UNCATEGORIZED_LABEL = "not_in_user_gene_list"

###########
# HELPERS #
###########


def _filter_empty_categories(
    gene_categories_df: pd.DataFrame, genes_in_rates: set[str], category_label: str
) -> pd.DataFrame:
    """Drop categories with no overlap with the rates matrix and log one warning each.

    Args:
        gene_categories_df: DataFrame indexed by gene ID with a ``category`` column.
        genes_in_rates: Set of gene IDs present in the expected mutations matrix.
        category_label: Human-readable name for the category type, used in warning messages
            (e.g. ``"user gene category"`` or ``"gnomAD LOEUF decile"``).

    Returns:
        Filtered copy of ``gene_categories_df`` with empty categories removed.
    """
    logger = logging.getLogger("logger")
    valid: list = []
    for category, cat_df in gene_categories_df.groupby("category"):
        if set(cat_df.index) & genes_in_rates:
            valid.append(category)
        else:
            logger.warning(f"No genes from the rates file belong to {category_label} '{category}'; skipping.")

    return gene_categories_df[gene_categories_df["category"].isin(valid)]


def _load_gene_infos(gff: Path) -> pd.DataFrame:
    """Extract chromosome and HGNC ID for every gene in a GFF database.

    Args:
        gff: Path to the GFF file (loaded via ``load_gff``).

    Returns:
        DataFrame indexed by gene ID with columns ``chrom`` and ``hgnc_id``.
        Genes without an HGNC identifier carry an empty string.
    """
    gff_db = load_gff(gff)

    chrom_map: dict[str, str] = {}
    hgnc_map: dict[str, str] = {}
    for gene in gff_db.features_of_type("gene"):
        chrom_map[gene.id] = gene.chrom
        attrs = dict(gene.attributes)
        hgnc_map[gene.id] = attrs["hgnc_id"][0] if "hgnc_id" in attrs else ""

    gene_infos = pd.DataFrame({"chrom": chrom_map, "hgnc_id": hgnc_map})
    return gene_infos


def _compute_observed_dnm_per_gene(dnm_snv_df: pd.DataFrame) -> pd.DataFrame:
    """Count observed de novo mutations per gene per consequence.

    Args:
        dnm_snv_df: SNV de novo mutation table with ``gene_id`` and
            ``consequence`` columns.

    Returns:
        DataFrame indexed by gene ID with one column per consequence and
        zero-filled missing entries.
    """
    observed: dict[str, dict[str, int]] = {}
    for gene_id, gene_df in dnm_snv_df.groupby("gene_id"):
        observed[gene_id] = {csq: csq_df.shape[0] for csq, csq_df in gene_df.groupby("consequence")}

    result = pd.DataFrame(observed).transpose()
    # reindex: adds consequence columns absent from every gene (fill_value=0)
    # fillna:  fills NaN for genes that had no mutations for a consequence that other genes did have
    result = result.reindex(columns=CONSEQUENCES, fill_value=0).fillna(0)
    return result


def _compute_expected_dnm_per_gene(
    rates_per_gene_per_csq_df: pd.DataFrame,
    gene_infos: pd.DataFrame,
    nmales: int,
    nfemales: int,
) -> pd.DataFrame:
    """Scale per-gene mutation rates to the cohort size.

    Args:
        rates_per_gene_per_csq_df: Sum of per-site mutation probabilities,
            indexed by gene ID, one column per consequence.
        gene_infos: DataFrame with ``chrom`` column, indexed by gene ID.
        nmales: Number of males in the cohort.
        nfemales: Number of females in the cohort.

    Returns:
        DataFrame of expected mutation counts, same shape as the input.
    """
    return rates_per_gene_per_csq_df.apply(
        lambda x: compute_expected_number_of_mutations_stats(x, gene_infos, nmales, nfemales),
        axis=1,
    ).copy()


############
# ANALYSIS #
############


def high_level_analysis(
    observed_dnm_per_gene_per_csq_df: pd.DataFrame,
    expected_mutations_per_gene_per_csq_df: pd.DataFrame,
    outdir: Path,
) -> None:
    """Compute genome-wide O/E ratios and write them to disk.

    Args:
        observed_dnm_per_gene_per_csq_df: Observed DNM counts per gene per consequence.
        expected_mutations_per_gene_per_csq_df: Expected DNM counts per gene per consequence.
        outdir: Directory where ``oe_ratios.json`` will be written.
    """

    logger = logging.getLogger("logger")
    logger.info(f"Computing observed/expected ratios for {expected_mutations_per_gene_per_csq_df.shape[0]} genes...")

    oe_ratios, oe_ratios_per_csq = compute_oe_ratios(
        observed_dnm_per_gene_per_csq_df, expected_mutations_per_gene_per_csq_df
    )
    _export_oe_ratios(oe_ratios, oe_ratios_per_csq, outdir)


def compute_oe_ratios(
    observed_dnm_per_gene_per_csq_df: pd.DataFrame,
    expected_mutations_per_gene_per_csq_df: pd.DataFrame,
) -> tuple[pd.Series, pd.DataFrame]:
    """Compute observed/expected ratios at gene level, globally and per consequence.

    Genes present in the expected table but absent from the observed table are
    treated as having zero observed mutations.

    Args:
        observed_dnm_per_gene_per_csq_df: Observed DNM counts per gene per consequence.
        expected_mutations_per_gene_per_csq_df: Expected DNM counts per gene per consequence.

    Returns:
        A pair ``(oe_ratios, oe_ratios_per_csq)`` where:

        - ``oe_ratios``: Series of aggregate O/E ratio per gene (all consequences summed).
        - ``oe_ratios_per_csq``: DataFrame of O/E ratio per gene per consequence.
    """
    # Pad observed with zeros for genes that had no mutations in the cohort
    genes_not_observed = set(expected_mutations_per_gene_per_csq_df.index) - set(observed_dnm_per_gene_per_csq_df.index)
    if genes_not_observed:
        padding = pd.DataFrame(0, index=list(genes_not_observed), columns=observed_dnm_per_gene_per_csq_df.columns)
        observed_dnm_per_gene_per_csq_df = pd.concat([observed_dnm_per_gene_per_csq_df, padding])

    oe_ratios = observed_dnm_per_gene_per_csq_df.sum(axis=1) / expected_mutations_per_gene_per_csq_df.sum(axis=1)
    oe_ratios_per_csq = observed_dnm_per_gene_per_csq_df / expected_mutations_per_gene_per_csq_df

    return oe_ratios, oe_ratios_per_csq


def category_analysis(
    observed_dnm_per_gene_per_csq_df: pd.DataFrame,
    expected_mutations_per_gene_per_csq_df: pd.DataFrame,
    gene_categories_df: pd.DataFrame,
    outdir: Path,
) -> None:
    """Plot and export O/E ratios per user-defined gene category for each consequence.

    Args:
        observed_dnm_per_gene_per_csq_df: Observed DNM counts per gene per consequence.
        expected_mutations_per_gene_per_csq_df: Expected DNM counts per gene per consequence.
        gene_categories_df: DataFrame indexed by gene ID with a ``category`` column.
        outdir: Directory where one PNG per consequence and a summary JSON will be written.
    """
    logger = logging.getLogger("logger")
    logger.info("Computing observed/expected ratios per user gene categories...")

    # Warn once per empty category before iterating over consequences
    genes_in_rates = set(expected_mutations_per_gene_per_csq_df.index)
    gene_categories_df = _filter_empty_categories(
        gene_categories_df, genes_in_rates, category_label="user gene category"
    )

    category_oe: dict[str, pd.DataFrame] = {}

    for csq in CONSEQUENCES:
        # Defensive: _compute_observed_dnm_per_gene guarantees all CONSEQUENCES columns,
        # but external callers may pass a narrower DataFrame.
        if csq not in observed_dnm_per_gene_per_csq_df.columns:
            logger.warning(f"Consequence '{csq}' not found in observed data; skipping category analysis.")
            continue

        oe_ratio_df = _oe_ratio_per_category(
            observed_dnm_per_gene_per_csq_df[[csq]],
            expected_mutations_per_gene_per_csq_df[[csq]],
            gene_categories_df,
        )
        category_oe[csq] = oe_ratio_df
        _plot_oe_ratio_per_gene_category(oe_ratio_df, csq, outdir / f"oe_ratio_per_category_{csq}.png")

    _export_category_oe_ratios(category_oe, outdir, "oe_ratios_per_category.json")


def _oe_ratio_per_category(
    df_observed: pd.DataFrame,
    df_expected: pd.DataFrame,
    df_gene_category: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate O/E ratios and 95 % Poisson confidence intervals by gene category.

    Args:
        df_observed: Observed DNM counts (one consequence column), indexed by gene ID.
        df_expected: Expected DNM counts (one consequence column), indexed by gene ID.
        df_gene_category: DataFrame indexed by gene ID with a ``category`` column.

    Returns:
        DataFrame with one row per category and columns ``oe_ratio``,
        ``nb_genes_in_category``, ``sum_observed``, ``sum_expected``,
        ``ci.lower``, ``ci.upper``, and ``category``.
    """
    s_expected = df_expected.sum(axis=1)
    s_observed = df_observed.sum(axis=1)
    observed_genes = set(s_observed.index)

    rows = []
    categorized_genes: set[str] = set()

    for category, category_df in df_gene_category.groupby("category"):
        genes = set(category_df.index) & set(s_expected.index)
        categorized_genes |= genes
        sum_observed = s_observed[list(genes & observed_genes)].sum()
        sum_expected = s_expected[list(genes)].sum()

        # Some genes may have zero expected mutations for a consequence (e.g. splicing in single-exon genes)
        if sum_expected == 0:
            continue

        rows.append(
            pd.Series(
                {
                    "oe_ratio": sum_observed / sum_expected,
                    "nb_genes_in_category": len(genes),
                    "sum_observed": sum_observed,
                    "sum_expected": sum_expected,
                },
                name=category,
            )
        )

    # Genes present in the rates matrix but absent from the category file
    uncategorized = set(s_expected.index) - categorized_genes
    if uncategorized:
        sum_observed = s_observed[list(uncategorized & set(s_observed.index))].sum()
        sum_expected = s_expected[list(uncategorized)].sum()
        if sum_expected > 0:
            rows.append(
                pd.Series(
                    {
                        "oe_ratio": sum_observed / sum_expected,
                        "nb_genes_in_category": len(uncategorized),
                        "sum_observed": sum_observed,
                        "sum_expected": sum_expected,
                    },
                    name=_UNCATEGORIZED_LABEL,
                )
            )

    df = pd.DataFrame(rows)

    # Poisson 95 % confidence intervals on the observed count, scaled to O/E space
    df["sum_observed"] = df["sum_observed"].astype(int)
    df["ci.lower"] = df.apply(lambda row: poisson.interval(0.95, row["sum_observed"])[0] / row["sum_expected"], axis=1)
    df["ci.upper"] = df.apply(lambda row: poisson.interval(0.95, row["sum_observed"])[1] / row["sum_expected"], axis=1)
    df["category"] = df.index

    return df


def _plot_oe_ratio_per_gene_category(
    df: pd.DataFrame,
    csq: str,
    output_path: Path,
    category_label: str = "gene category",
) -> None:
    """Save a scatter plot with error bars of O/E ratio per category.

    Args:
        df: Output of :func:`_oe_ratio_per_category`.
        csq: Consequence label used in the plot title.
        output_path: Destination PNG file.
        category_label: Human-readable name for the category axis, used in the
            x-axis label and plot title (e.g. ``"gnomAD missense z-score decile"``).
    """
    # User-defined categories sorted alphabetically; uncategorized always last
    user_cats = sorted(c for c in df["category"] if c != _UNCATEGORIZED_LABEL)
    categories = user_cats + ([_UNCATEGORIZED_LABEL] if _UNCATEGORIZED_LABEL in df["category"].values else [])

    df = df.set_index("category").loc[categories].reset_index()
    x_pos = np.arange(len(categories))

    colors = ["lightgray" if c == _UNCATEGORIZED_LABEL else "steelblue" for c in categories]
    ys = df["oe_ratio"].values
    yerr = [ys - df["ci.lower"].values, df["ci.upper"].values - ys]

    fig, ax = plt.subplots(figsize=(max(6, len(categories) * 1.4), 6))

    ax.scatter(x_pos, ys, s=80, color=colors, zorder=3)
    ax.errorbar(x_pos, ys, yerr=yerr, fmt="none", ecolor="black", capsize=4, elinewidth=1, zorder=2)
    ax.axhline(1, color="red", linestyle="--", linewidth=1, zorder=1)

    # Gene count on a second line below each category label
    tick_labels = [f"{c}\n(n={int(n)})" for c, n in zip(categories, df["nb_genes_in_category"])]
    ax.set_xticks(x_pos)
    ax.set_xticklabels(tick_labels, rotation=30, ha="right")
    ax.set_xlabel(category_label)
    ax.set_ylabel("O/E ratio")
    ax.set_title(f"Observed/expected {csq} DNM ratio per {category_label}")
    ax.set_ylim(bottom=0)
    ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def gnomad_constraints_analysis(
    observed_dnm_per_gene_per_csq_df: pd.DataFrame,
    expected_mutations_per_gene_per_csq_df: pd.DataFrame,
    gnomad_constraints_df: pd.DataFrame,
    outdir: Path,
) -> None:
    """Plot and export O/E ratios stratified by gnomAD constraint metrics.

    Genes are stratified into deciles of:

    - missense z-score (higher = more missense-constrained)
    - LOEUF (lower = more LoF-constrained)

    Gene identifiers in the gnomAD file are matched to the versioned IDs used in the
    rates matrix via :func:`~denovowest.utils.helpers.ensembl_ids_map_version`.
    Genes that cannot be matched are logged as a warning and excluded.

    Args:
        observed_dnm_per_gene_per_csq_df: Observed DNM counts per gene per consequence.
        expected_mutations_per_gene_per_csq_df: Expected DNM counts per gene per consequence.
        gnomad_constraints_df: gnomAD constraint table with at least ``gene_id``,
            ``mis.z_score``, and ``lof.oe_ci.upper`` columns, indexed arbitrarily.
        outdir: Directory where plots and JSON files will be written.
    """
    logger = logging.getLogger("logger")
    logger.info("Computing observed/expected ratios per gnomad constraints bins...")

    # TODO : should be better handled to look for transcripts provided in gff
    # Keep only genes with ENSG IDs and a unique entry in the gnomAD constraints file, preferring MANE Select or canonical transcripts where available
    gnomad_constraints_df = gnomad_constraints_df.loc[gnomad_constraints_df["gene_id"].str.startswith("ENSG")].copy()
    gnomad_constraints_df = gnomad_constraints_df.loc[
        gnomad_constraints_df.mane_select | gnomad_constraints_df.canonical
    ].copy()
    gnomad_constraints_df = gnomad_constraints_df.loc[~gnomad_constraints_df.gene_id.duplicated()]

    # Remap gnomAD gene IDs to match the versioned IDs used in the rates matrix
    id_mapping = ensembl_ids_map_version(
        list(gnomad_constraints_df.gene_id),
        list(expected_mutations_per_gene_per_csq_df.index),
    )
    gnomad_constraints_df.index = gnomad_constraints_df.gene_id.map(id_mapping)
    gnomad_constraints_df = gnomad_constraints_df.loc[
        gnomad_constraints_df.index.isin(expected_mutations_per_gene_per_csq_df.index)
    ].copy()

    # Log genes present in the rates matrix but absent from gnomAD
    n_unmatched = len(set(expected_mutations_per_gene_per_csq_df.index) - set(gnomad_constraints_df.index))
    if n_unmatched:
        logger.warning(
            f"{n_unmatched} gene(s) from the rates matrix have no matching entry in the "
            f"gnomAD constraints file and will be excluded from the stratified analysis."
        )

    # Add decile bins for missense z-score and LOEUF
    gnomad_constraints_df["mis.z_score.decile"] = (
        pd.qcut(gnomad_constraints_df["mis.z_score"], 10, labels=False) + 1
    ).astype(int)
    gnomad_constraints_df["lof.oe_ci.upper.decile"] = (
        pd.qcut(gnomad_constraints_df["lof.oe_ci.upper"], 10, labels=False) + 1
    ).astype(int)

    # Build and filter category DataFrames once — they are gene-level and consequence-independent
    genes_in_rates = set(expected_mutations_per_gene_per_csq_df.index)
    mis_cat_df = _filter_empty_categories(
        gnomad_constraints_df[["mis.z_score.decile"]].rename(columns={"mis.z_score.decile": "category"}),
        genes_in_rates,
        category_label="gnomAD missense z-score decile",
    )
    loeuf_cat_df = _filter_empty_categories(
        gnomad_constraints_df[["lof.oe_ci.upper.decile"]].rename(columns={"lof.oe_ci.upper.decile": "category"}),
        genes_in_rates,
        category_label="gnomAD LOEUF decile",
    )

    missense_oe: dict[str, pd.DataFrame] = {}
    loeuf_oe: dict[str, pd.DataFrame] = {}

    for csq in CONSEQUENCES:
        if csq not in observed_dnm_per_gene_per_csq_df.columns:
            logger.warning(f"Consequence '{csq}' not found in observed data; skipping gnomAD constraints analysis.")
            continue

        oe_mis = _oe_ratio_per_category(
            observed_dnm_per_gene_per_csq_df[[csq]], expected_mutations_per_gene_per_csq_df[[csq]], mis_cat_df
        )
        missense_oe[csq] = oe_mis
        _plot_oe_ratio_per_gene_category(
            oe_mis,
            csq,
            outdir / f"oe_ratio_per_gnomad_missense_bins_{csq}.png",
            category_label="gnomAD missense z-score decile",
        )

        oe_loeuf = _oe_ratio_per_category(
            observed_dnm_per_gene_per_csq_df[[csq]], expected_mutations_per_gene_per_csq_df[[csq]], loeuf_cat_df
        )
        loeuf_oe[csq] = oe_loeuf
        _plot_oe_ratio_per_gene_category(
            oe_loeuf,
            csq,
            outdir / f"oe_ratio_per_gnomad_loeuf_bins_{csq}.png",
            category_label="gnomAD LOEUF decile",
        )

    _export_category_oe_ratios(missense_oe, outdir, "oe_ratios_per_gnomad_missense_bins.json")
    _export_category_oe_ratios(loeuf_oe, outdir, "oe_ratios_per_gnomad_loeuf_bins.json")


##########
# OUTPUT #
##########


def _export_category_oe_ratios(
    oe_ratios_by_csq: dict[str, pd.DataFrame],
    outdir: Path,
    filename: str,
) -> None:
    """Serialise per-category O/E ratios to JSON, keyed by consequence then category.

    Args:
        oe_ratios_by_csq: Mapping from consequence name to the DataFrame returned by
            :func:`_oe_ratio_per_category`.
        outdir: Directory where the file will be written.
        filename: Output file name (e.g. ``oe_ratios_per_category.json``).
    """
    export_cols = ["oe_ratio", "nb_genes_in_category", "sum_observed", "sum_expected", "ci.lower", "ci.upper"]
    result = {
        csq: df.drop(columns="category")[export_cols].to_dict(orient="index") for csq, df in oe_ratios_by_csq.items()
    }
    with open(outdir / filename, "w") as f:
        json.dump(result, f, indent=2)


def _export_oe_ratios(oe_ratios: pd.Series, oe_ratios_per_csq: pd.DataFrame, outdir: Path) -> None:
    """Serialise O/E ratios to JSON, merging per-consequence and aggregate values.

    The output JSON maps each gene ID to a dict whose keys are consequence names
    plus ``"all"`` for the aggregate ratio.

    Args:
        oe_ratios: Aggregate O/E ratio per gene (all consequences combined).
        oe_ratios_per_csq: Per-consequence O/E ratio per gene.
        outdir: Directory where ``oe_ratios.json`` will be written.
    """
    per_csq_dict = oe_ratios_per_csq.to_dict(orient="index")
    aggregate_dict = {gene_id: {"all": ratio} for gene_id, ratio in oe_ratios.items()}

    merged: dict[str, dict] = {}
    for key in per_csq_dict.keys() | aggregate_dict.keys():
        merged[key] = {**per_csq_dict.get(key, {}), **aggregate_dict.get(key, {})}

    with open(outdir / "oe_ratios.json", "w") as f:
        json.dump(merged, f, indent=2)


#######
# CLI #
#######


@click.command()
@click.argument("dnm", type=click.Path(exists=True, path_type=Path))
@click.argument("rates_per_gene_per_csq", type=click.Path(exists=True, path_type=Path))
@click.argument("gff", type=click.Path(exists=True, path_type=Path))
@click.argument("nmales", type=int)
@click.argument("nfemales", type=int)
@click.argument("outdir", type=click.Path(exists=True, path_type=Path))
@click.option("--gnomad-constraints", type=click.Path(path_type=Path))
@click.option("--gene-categories", type=click.Path(path_type=Path))
def main(
    dnm: Path,
    rates_per_gene_per_csq: Path,
    gff: Path,
    nmales: int,
    nfemales: int,
    outdir: Path,
    gnomad_constraints: Optional[Path],
    gene_categories: Optional[Path],
) -> None:
    """Compute observed/expected DNM ratios per gene and per consequence.

    \b
    Arguments:
        DNM:                    Tab-delimited de novo mutation file.
        RATES_PER_GENE_PER_CSQ: Per-gene per-consequence mutation rate matrix (TSV).
        GFF:                    Genome annotation file used to extract gene metadata.
        NMALES:                 Number of males in the cohort.
        NFEMALES:               Number of females in the cohort.
        OUTDIR:                 Directory where output files will be written.
    """
    init_log()

    # Read DNM file; keep only SNVs with a consequence of interest
    dnm_df = pd.read_csv(dnm, sep="\t")
    dnm_snv_df = dnm_df.loc[(dnm_df["ref"].str.len() == 1) & (dnm_df["alt"].str.len() == 1)]
    dnm_snv_df = dnm_snv_df.loc[dnm_snv_df["consequence"].isin(CONSEQUENCES)]

    # Load the per-gene per-consequence mutation rate matrix
    rates_per_gene_per_csq_df = pd.read_csv(rates_per_gene_per_csq, sep="\t", index_col=0)
    rates_per_gene_per_csq_df = rates_per_gene_per_csq_df.loc[:, CONSEQUENCES]

    # Extract gene chromosome and HGNC metadata from the GFF
    gene_infos = _load_gene_infos(gff)

    # Count observed DNMs and compute expected counts scaled to the cohort
    observed_dnm_per_gene_per_csq_df = _compute_observed_dnm_per_gene(dnm_snv_df)
    expected_mutations_per_gene_per_csq_df = _compute_expected_dnm_per_gene(
        rates_per_gene_per_csq_df, gene_infos, nmales, nfemales
    )

    # Get O/E ratios genome wide
    high_level_analysis(observed_dnm_per_gene_per_csq_df, expected_mutations_per_gene_per_csq_df, outdir)

    # Get O/E ratios per gene categories
    if gene_categories:
        gene_categories_df = pd.read_csv(gene_categories, sep="\t", index_col=0)
        category_analysis(
            observed_dnm_per_gene_per_csq_df,
            expected_mutations_per_gene_per_csq_df,
            gene_categories_df,
            outdir,
        )

    # Get O/E ratios per gnomad constraint bins
    if gnomad_constraints:
        gnomad_constraints_df = pd.read_csv(gnomad_constraints, sep="\t", index_col=0)
        gnomad_constraints_analysis(
            observed_dnm_per_gene_per_csq_df,
            expected_mutations_per_gene_per_csq_df,
            gnomad_constraints_df,
            outdir,
        )


if __name__ == "__main__":
    main()
