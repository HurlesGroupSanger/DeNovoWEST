#!/usr/bin/env python

import json
from pathlib import Path
from typing import Optional

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import poisson

from denovowest.utils.helpers import compute_expected_number_of_mutations_stats
from denovowest.utils.io_helpers import load_gff
from denovowest.utils.params import ConsequenceGroups

##############
# CONSTANTS  #
##############

CONSEQUENCES = ConsequenceGroups.OE_RATIO

###########
# HELPERS #
###########


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
    # Reindex to ensure all tracked consequences are always present as columns,
    # even when none of the input DNMs carry that consequence
    result = result.reindex(columns=CONSEQUENCES, fill_value=0)
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
    gene_infos: pd.DataFrame,
    gene_categories_df: pd.DataFrame,
    outdir: Path,
) -> None:
    """Plot O/E ratios per gene category for each consequence.

    Args:
        observed_dnm_per_gene_per_csq_df: Observed DNM counts per gene per consequence.
        expected_mutations_per_gene_per_csq_df: Expected DNM counts per gene per consequence.
        gene_infos: Gene metadata (not used directly, reserved for future use).
        gene_categories_df: DataFrame indexed by gene ID with a ``category`` column.
        outdir: Directory where one PNG per consequence will be written.
    """
    for csq in CONSEQUENCES:
        try:
            observed_dnm_per_gene_df = observed_dnm_per_gene_per_csq_df[[csq]]
            expected_mutations_per_gene_df = expected_mutations_per_gene_per_csq_df[[csq]]
        except KeyError:
            print(f"Consequence {csq} not found in the data, skipping category analysis for this consequence.")
            continue

        oe_ratio_df = _oe_ratio_per_category(
            observed_dnm_per_gene_df, expected_mutations_per_gene_df, gene_categories_df
        )
        _plot_oe_ratio_per_gene_category(oe_ratio_df, csq, outdir / f"oe_ratio_per_category_{csq}.png")


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

    rows = []
    for category, min_df in df_gene_category.groupby("category"):
        genes = set(min_df.index) & set(s_expected.index)
        sum_observed = s_observed[list(genes & set(s_observed.index))].sum()
        sum_expected = s_expected[list(genes)].sum()

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

    df = pd.DataFrame(rows)

    # Poisson 95 % confidence intervals on the observed count, scaled to O/E space
    df["sum_observed"] = df["sum_observed"].astype(int)
    df["ci.lower"] = df.apply(lambda row: poisson.interval(0.95, row["sum_observed"])[0] / row["sum_expected"], axis=1)
    df["ci.upper"] = df.apply(lambda row: poisson.interval(0.95, row["sum_observed"])[1] / row["sum_expected"], axis=1)
    df["category"] = df.index

    return df


def _plot_oe_ratio_per_gene_category(df: pd.DataFrame, csq: str, output_path: Path) -> None:
    """Save a scatter plot with error bars of O/E ratio per gene category.

    Args:
        df: Output of :func:`_oe_ratio_per_category`.
        csq: Consequence label used in the plot title.
        output_path: Destination PNG file.
    """
    categories = df["category"].unique()
    x_pos = np.arange(len(categories))

    ys = df["oe_ratio"].values
    yerr = [ys - df["ci.lower"].values, df["ci.upper"].values - ys]

    plt.figure(figsize=(14, 8))
    plt.scatter(x_pos, ys, s=80, color="steelblue")
    plt.errorbar(x_pos, ys, yerr=yerr, fmt="none", ecolor="black", capsize=4, elinewidth=1)
    plt.axhline(1, color="red", linestyle="--")
    plt.xticks(x_pos, categories)
    plt.xlabel("Gene category")
    plt.ylabel("O/e ratio")
    plt.title(f"Observed/expected {csq} DNM ratio per gene category")
    plt.ylim(bottom=0)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def gnomad_constraints_analysis(
    observed_dnm_per_gene_per_csq_df: pd.DataFrame,
    expected_mutations_per_gene_per_csq_df: pd.DataFrame,
    gnomad_constraints_df: pd.DataFrame,
    outdir: Path,
) -> None:
    # TODO: implement gnomad constraints analysis
    pass


##########
# OUTPUT #
##########


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

    high_level_analysis(observed_dnm_per_gene_per_csq_df, expected_mutations_per_gene_per_csq_df, outdir)

    if gene_categories:
        gene_categories_df = pd.read_csv(gene_categories, sep="\t", index_col=0)
        category_analysis(
            observed_dnm_per_gene_per_csq_df,
            expected_mutations_per_gene_per_csq_df,
            gene_infos,
            gene_categories_df,
            outdir,
        )

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
