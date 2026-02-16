#!/usr/bin/env python

import click
import pandas as pd
from scipy.stats import combine_pvalues

# dne_file_all = (
#     "../input/merged_all_dne_test_ppv_2020_03_09.tab"
# )
# dne_file_mis = (
#     "../input/merged_mis_dne_test_ppv_2020_03_09.tab"
# )
# dnn_file = "../input/denovonear_out_missense_31058_ntrios_2019_05_15.txt"
# out_file = "out.tsv"


@click.command()
@click.argument("dne_file_all", type=click.Path(exists=True))
@click.argument("dne_file_mis", type=click.Path(exists=True))
@click.argument("dnn_file", type=click.Path(exists=True))
@click.argument("out_file", type=click.Path())
def combine_dne_dnn(dne_file_all, dne_file_mis, dnn_file, out_file):
    """
    Combine p-values from different sources using Fisher's method (log-likelihood ratio test).
    Extract the minimum of the combined p-values.

    Args:
        dne_file_all (str): Path to the file containing DNW results for all variants.
        dne_file_mis (str): Path to the file containing DNW results for missense variants.
        dnn_file (str): Path to the file containing DNN results.
        out_file (str): Path to the output file.

    Returns:
        None
    """

    # Load DNW results
    dne_all = pd.read_csv(dne_file_all, sep="\t")
    dne_all.columns = ["gene_id", "expected_all", "observed_all", "enrichment_all_pval"]
    dne_mis = pd.read_csv(dne_file_mis, sep="\t")
    dne_mis.columns = [
        "gene_id",
        "expected_mis",
        "observed_mis",
        "enrichment_mis_pval",
    ]

    # Load DNN results
    dnn = pd.read_csv(dnn_file, sep="\t")
    dnn = dnn[dnn["mutation_category"] == "missense"]
    dnn.rename(
        {"probability": "clustering_pval", "events_n": "nb_mis_variants", "dist": "clustering_dist"},
        axis=1,
        inplace=True,
    )

    # Merge dataframes
    merged_data = pd.merge(dne_all, dnn, on="gene_id", how="left")
    merged_data = pd.merge(merged_data, dne_mis, on="gene_id", how="outer")

    # Combine p-values using Fisher's method (log-likelihood ratio test)
    merged_data["combined_mis_pval"] = merged_data[["clustering_pval", "enrichment_mis_pval"]].apply(
        lambda pvals: combine_pvalues(pvals, method="fisher")[1], axis=1
    )

    # Take minimum of combined and enrichment p-values
    merged_data["min_pval"] = merged_data[["combined_mis_pval", "enrichment_all_pval"]].apply(
        lambda pvals: min(pvals.dropna(), default=1), axis=1
    )

    # Write to file
    merged_data = merged_data[
        [
            "gene_id",
            "min_pval",
            "combined_mis_pval",
            "enrichment_all_pval",
            "enrichment_mis_pval",
            "clustering_pval",
            "observed_all",
            "expected_all",
            "observed_mis",
            "expected_mis",
            "nb_mis_variants",
            "clustering_dist",
        ]
    ]
    merged_data.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    combine_dne_dnn()
