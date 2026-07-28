#!/usr/bin/env python

import click
import pandas as pd
from scipy.stats import combine_pvalues


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

    # Load DNW results — rename from simulation output column names to analysis-specific names
    dne_all = pd.read_csv(dne_file_all, sep="\t").rename(
        columns={
            "expected_score": "expected_all",
            "observed_score": "observed_all",
            "p-value": "enrichment_all_pval",
        }
    )
    dne_mis = pd.read_csv(dne_file_mis, sep="\t").rename(
        columns={
            "expected_score": "expected_mis",
            "observed_score": "observed_mis",
            "p-value": "enrichment_mis_pval",
        }
    )

    # Extract gene_symbol from dne_all before the merge — it is the authoritative source because
    # dne_mis only covers genes with missense mutations and may not include all genes
    gene_symbol_df = None
    if "gene_symbol" in dne_all.columns:
        gene_symbol_df = dne_all[["gene_id", "gene_symbol"]].drop_duplicates("gene_id")
    dne_all = dne_all.drop(columns=["gene_symbol"], errors="ignore")
    dne_mis = dne_mis.drop(columns=["gene_symbol"], errors="ignore")

    # Load DNN results
    dnn = pd.read_csv(dnn_file, sep="\t")
    dnn = dnn[dnn["mutation_category"] == "missense"]
    dnn.rename(
        {"probability": "clustering_pval", "events_n": "clustering_nb_mis_variants", "dist": "clustering_dist"},
        axis=1,
        inplace=True,
    )
    # If the user combined 3D and linear clustering we want to keep track of the mode for each gene
    # — must check before subsetting columns or the column will be dropped
    if "mode" in dnn.columns:
        dnn.rename({"mode": "clustering_mode"}, axis=1, inplace=True)

    cols = ["gene_id", "clustering_pval", "clustering_nb_mis_variants", "clustering_dist"]
    if "clustering_mode" in dnn.columns:
        cols.append("clustering_mode")
    dnn = dnn[cols]

    # Merge dataframes
    merged_data = pd.merge(dne_all, dnn, on="gene_id", how="left")
    merged_data = pd.merge(merged_data, dne_mis, on="gene_id", how="outer")

    # Re-attach gene_symbol from the dne_all lookup (left join so genes absent from dne_all get NaN)
    if gene_symbol_df is not None:
        merged_data = merged_data.merge(gene_symbol_df, on="gene_id", how="left")

    # Combine p-values using Fisher's method (log-likelihood ratio test)
    merged_data["combined_mis_pval"] = merged_data[["clustering_pval", "enrichment_mis_pval"]].apply(
        lambda pvals: combine_pvalues(pvals, method="fisher")[1], axis=1
    )

    # Take minimum of combined and enrichment p-values
    merged_data["min_pval"] = merged_data[["combined_mis_pval", "enrichment_all_pval"]].apply(
        lambda pvals: min(pvals.dropna(), default=1), axis=1
    )

    # Write to file
    columns_to_keep = ["gene_id"]
    if "gene_symbol" in merged_data.columns:
        columns_to_keep.append("gene_symbol")
    columns_to_keep += [
        "min_pval",
        "enrichment_all_pval",
        "combined_mis_pval",
        "enrichment_mis_pval",
        "clustering_pval",
        "observed_all",
        "expected_all",
        "observed_mis",
        "expected_mis",
        "clustering_nb_mis_variants",
        "clustering_dist",
    ]
    if "clustering_mode" in merged_data.columns:
        columns_to_keep.append("clustering_mode")

    merged_data = merged_data[columns_to_keep]
    merged_data.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    combine_dne_dnn()
