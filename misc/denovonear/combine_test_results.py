#!/usr/bin/env python
import click
import pandas as pd
import numpy as np


@click.command()
@click.argument("clustering_linear", type=click.Path(exists=True))
@click.argument("clustering_3d", type=click.Path(exists=True))
@click.argument("out", type=str)
def combine_test_results(clustering_linear, clustering_3d, out):
    """
    Combine linear and 3D clustering test.
    When a 3D results is available we pick it, otherwise we pick the linear results.

    Args:
        clustering_linear (str): linear clustering results
        clustering_3d (str): 3d clustering results
        out (str): combined clustering results
    """
    clustering_linear_df = pd.read_csv(clustering_linear, sep="\t", index_col=0)
    clustering_linear_df = clustering_linear_df.loc[clustering_linear_df.mutation_category == "missense"]

    clustering_3d_df = pd.read_csv(clustering_3d, sep="\t", index_col=0)
    clustering_3d_df = clustering_3d_df.loc[clustering_3d_df.mutation_category == "missense"]

    merged_clustering = clustering_linear_df.merge(
        clustering_3d_df, left_index=True, right_index=True, how="outer", suffixes=["_linear", "_3d"]
    )

    list_series = list()
    for gene_id, row in merged_clustering.iterrows():
        if np.isnan(row.probability_3d):
            symbol = row.gene_symbol_linear
            events_n = row.events_n_linear
            dist = row.dist_linear
            probability = row.probability_linear
        else:
            symbol = row.gene_symbol_3d
            events_n = row.events_n_3d
            dist = row.dist_3d
            probability = row.probability_3d

        s = pd.Series(
            {
                "gene_symbol": symbol,
                "mutation_category": "missense",
                "events_n": events_n,
                "dist": dist,
                "probability": probability,
            }
        )
        s.name = gene_id
        list_series.append(s)

    res_df = pd.DataFrame(list_series)
    res_df.index.name = "gene_id"

    res_df.to_csv(out, sep="\t")


if __name__ == "__main__":
    combine_test_results()
