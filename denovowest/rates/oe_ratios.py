import json
from pathlib import Path

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import poisson

from denovowest.utils.helpers import compute_expected_number_of_mutations_stats
from denovowest.utils.io_helpers import load_gff

CONSEQUENCES = [
    "synonymous",
    "missense",
    "stop_gained",
    "stop_retained",
    "splice_region",
    "splice_donor",
    "splice_acceptor",
    "start_lost",
    "stop_lost",
]

###########
# HELPERS #
###########


def get_gene_infos(gff):

    gff_db = load_gff(gff)

    gene_chrom_dict = dict()
    gene_hgnc_dict = dict()
    for gene in gff_db.features_of_type("gene"):
        gene_chrom_dict[gene.id] = gene.chrom
        if "hgnc_id" not in dict(gene.attributes).keys():
            gene_hgnc_dict[gene.id] = ""
        else:
            gene_hgnc_dict[gene.id] = dict(gene.attributes)["hgnc_id"][0]

    gene_infos = pd.DataFrame([gene_chrom_dict, gene_hgnc_dict]).transpose()
    gene_infos.columns = ["chrom", "hgnc_id"]

    return gene_infos


def get_observed_number_of_mutations_per_gene_per_consequence(dnm_snv_df):

    observed_dnm_per_gene_per_csq = dict()
    for gene_id, min_df in dnm_snv_df.groupby("gene_id"):
        observed_dnm_per_gene_per_csq[gene_id] = dict()
        for consequence, min_csq_df in min_df.groupby("consequence"):
            observed_dnm_per_gene_per_csq[gene_id][consequence] = min_csq_df.shape[0]

    observed_dnm_per_gene_per_csq_df = pd.DataFrame(observed_dnm_per_gene_per_csq).transpose()
    observed_dnm_per_gene_per_csq_df.fillna(0, inplace=True)

    return observed_dnm_per_gene_per_csq_df


def get_expected_number_of_mutations_per_gene_per_consequence(rates_per_gene_per_csq_df, gene_infos, nmales, nfemales):

    expected_mutations_per_gene_per_csq_df = rates_per_gene_per_csq_df.apply(
        lambda x: compute_expected_number_of_mutations_stats(x, gene_infos, nmales, nfemales), axis=1
    ).copy()

    return expected_mutations_per_gene_per_csq_df


############
# ANALYSIS #
############


def high_level_analysis(observed_dnm_per_gene_per_csq_df, expected_mutations_per_gene_per_csq_df, outdir):

    # Compute O/E ratios
    oe_ratios, oe_ratios_per_csq = compute_oe_ratios(
        observed_dnm_per_gene_per_csq_df, expected_mutations_per_gene_per_csq_df
    )

    export_oe_ratios(oe_ratios, oe_ratios_per_csq, outdir)


def compute_oe_ratios(observed_dnm_per_gene_per_csq_df, expected_mutations_per_gene_per_csq_df):

    genes_not_observed = set(expected_mutations_per_gene_per_csq_df.index) - set(observed_dnm_per_gene_per_csq_df.index)

    observed_dnm_per_gene_per_csq_df = pd.concat(
        [
            observed_dnm_per_gene_per_csq_df,
            pd.DataFrame(0, index=list(genes_not_observed), columns=observed_dnm_per_gene_per_csq_df.columns),
        ]
    )

    oe_ratios = observed_dnm_per_gene_per_csq_df.sum(axis=1) / expected_mutations_per_gene_per_csq_df.sum(axis=1)
    oe_ratios_per_csq = observed_dnm_per_gene_per_csq_df / expected_mutations_per_gene_per_csq_df

    return oe_ratios, oe_ratios_per_csq


def category_analysis(
    observed_dnm_per_gene_per_csq_df, expected_mutations_per_gene_per_csq_df, gene_infos, gene_categories_df, outdir
):

    for csq in CONSEQUENCES:

        try:
            observed_dnm_per_gene_df = observed_dnm_per_gene_per_csq_df[[csq]]
            expected_mutations_per_gene_df = expected_mutations_per_gene_per_csq_df[[csq]]
        except KeyError:
            print(f"Consequence {csq} not found in the data, skipping category analysis for this consequence.")
            continue

        oe_ratio_df = oe_ratio_per_category(
            observed_dnm_per_gene_df, expected_mutations_per_gene_df, gene_categories_df
        )

        plot_oe_ratio_per_gene_category(oe_ratio_df, csq, outdir / f"oe_ratio_per_category_{csq}.png")


def oe_ratio_per_category(df_observed, df_expected, df_gene_category):

    s_expected = df_expected.sum(axis=1)
    s_observed = df_observed.sum(axis=1)

    list_s = list()
    for category, min_df_gene_category in df_gene_category.groupby("category"):

        genes_in_category = set(min_df_gene_category.index) & set(s_expected.index)
        nb_genes_in_category = len(genes_in_category)
        sum_observed = s_observed[list(set(s_observed.index) & set(genes_in_category))].sum()
        sum_expected = s_expected[list(genes_in_category)].sum()
        oe_ratio = sum_observed / sum_expected

        s = pd.Series(
            {
                "oe_ratio": oe_ratio,
                "nb_genes_in_category": nb_genes_in_category,
                "sum_observed": sum_observed,
                "sum_expected": sum_expected,
            },
            name=category,
        )
        list_s.append(s)

    df = pd.DataFrame(list_s)

    # Compute Confidence interval
    df.sum_observed = df.sum_observed.astype(int)
    df["ci.lower"] = df.apply(lambda row: poisson.interval(0.95, row["sum_observed"])[0] / row["sum_expected"], axis=1)
    df["ci.upper"] = df.apply(lambda row: poisson.interval(0.95, row["sum_observed"])[1] / row["sum_expected"], axis=1)

    df["category"] = df.index

    return df


def plot_oe_ratio_per_gene_category(df, csq, output_path):
    categories = df["category"].unique()
    x_base = np.arange(len(categories))

    plt.figure(figsize=(14, 8))

    ys = df["oe_ratio"].values
    yerr = [ys - df["ci.lower"].values, df["ci.upper"].values - ys]
    plt.scatter(x_base, ys, s=80, color="steelblue")
    plt.errorbar(x_base, ys, yerr=yerr, fmt="none", ecolor="black", capsize=4, elinewidth=1)

    plt.axhline(1, color="red", linestyle="--")
    plt.xticks(x_base, categories)
    plt.xlabel("Gene category")
    plt.ylabel("O/e ratio")
    plt.title(f"Observed/expected {csq} DNM ratio per gene category")
    plt.ylim(bottom=0)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def gnomad_constraints_analysis(
    observed_dnm_per_gene_per_csq_df, expected_mutations_per_gene_per_csq_df, gnomad_constraints_df, outdir
):
    # TODO: implement gnomad constraints analysis
    pass


##########
# OUTPUT #
##########


def export_oe_ratios(oe_ratios, oe_ratios_per_csq, outdir):

    oe_ratios_dict = oe_ratios_per_csq.to_dict(orient="index")
    oe_ratios_per_csq_dict = {gene_id: {"all": oe_ratio} for (gene_id, oe_ratio) in dict(oe_ratios).items()}

    merged_oe_ratios = {}
    for key in oe_ratios_dict.keys() | oe_ratios_per_csq_dict.keys():
        merged_oe_ratios[key] = {**oe_ratios_dict.get(key, {}), **oe_ratios_per_csq_dict.get(key, {})}

    with open(outdir / "oe_ratios.json", "w") as f:
        json.dump(merged_oe_ratios, f, indent=2)


@click.command()
@click.argument("dnm", type=click.Path(exists=True, path_type=Path))
@click.argument("rates_per_gene_per_csq", type=click.Path(exists=True, path_type=Path))
@click.argument("gff", type=click.Path(exists=True, path_type=Path))
@click.argument("nmales", type=int)
@click.argument("nfemales", type=int)
@click.argument("outdir", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--gnomad-constraints",
    type=click.Path(path_type=Path),
)
@click.option(
    "--gene-categories",
    type=click.Path(path_type=Path),
)
def main(
    dnm,
    rates_per_gene_per_csq,
    gff,
    nmales,
    nfemales,
    outdir,
    gnomad_constraints,
    gene_categories,
):
    """ """

    # Read DNM file and keep only SNVs with a consequence of interest
    dnm_df = pd.read_csv(dnm, sep="\t")
    dnm_snv_df = dnm_df.loc[
        dnm_df.apply(lambda row: True if (len(row.ref) == len(row.alt)) and (len(row.ref) == 1) else False, axis=1)
    ]
    dnm_snv_df = dnm_snv_df.loc[dnm_snv_df.consequence.isin(CONSEQUENCES)]

    # Read the rates per gene per consequence
    rates_per_gene_per_csq_df = pd.read_csv(rates_per_gene_per_csq, sep="\t", index_col=0)
    rates_per_gene_per_csq_df = rates_per_gene_per_csq_df.loc[:, CONSEQUENCES]

    # Read the GFF to extract gene information
    gene_infos = get_gene_infos(gff)

    # Compute the observed number of mutations per gene per consequence for the cohort
    observed_dnm_per_gene_per_csq_df = get_observed_number_of_mutations_per_gene_per_consequence(dnm_snv_df)

    # Compute the expected number of mutations per gene per consequence for the cohort (rates are not adjusted to the cohort size yet)
    expected_mutations_per_gene_per_csq_df = get_expected_number_of_mutations_per_gene_per_consequence(
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
            observed_dnm_per_gene_per_csq_df, expected_mutations_per_gene_per_csq_df, gnomad_constraints_df, outdir
        )


if __name__ == "__main__":
    main()
