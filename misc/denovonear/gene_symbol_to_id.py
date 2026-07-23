#!/usr/bin/env python
import click
import gffutils
import pandas as pd


@click.command()
@click.argument("dnn_results", type=click.Path(exists=True))
@click.argument("gff", type=click.Path(exists=True))
@click.argument("out", type=str)
def gene_symbol_to_id(dnn_results, gff, out):
    """
    Add ENSEMBL identifier to DNN results so we can merge it with DNW results.

    Args:
        dnn_results (_type_): _description_
        gff (_type_): _description_
        out (_type_): _description_
    """

    # Load DNM file
    dnn_df = pd.read_csv(dnn_results, sep="\t")
    dnn_df.rename({"gene_id": "gene_symbol"}, axis=1, inplace=True)

    # Load GFF
    if not gff.endswith("db"):
        gff_db = gffutils.create_db(gff, f"{gff}.db", merge_strategy="create_unique")
    else:
        gff_db = gffutils.FeatureDB(gff)

    # DeNovoNear uses gene name rather than gene identifiers
    gene_dict = dict()
    gene_duplicate_symbol = list()
    for gene in gff_db.features_of_type("gene", order_by="start"):
        gene_name = gene["Name"][0]
        gene_id = gene["ID"][0]
        if gene_name not in gene_dict:
            gene_dict[gene_name] = gene_id
        elif gene_name in list(dnn_df.gene_symbol):
            gene_duplicate_symbol.append((gene_id, gene_name))

    # Get the gene identifier
    dnn_df["gene_id"] = dnn_df["gene_symbol"].map(gene_dict)

    # Export Denovonear with ensembl IDs
    dnn_df = dnn_df[["gene_id", "gene_symbol", "mutation_category", "events_n", "dist", "probability"]]

    # Handle duplicated gene symbols
    for gene_id, gene_symbol in gene_duplicate_symbol:
        s = dnn_df.loc[(dnn_df["gene_symbol"] == gene_symbol) & (dnn_df.mutation_category == "missense")].iloc[0]
        new_s = s.copy()
        new_s["gene_id"] = gene_id
        dnn_df.loc[len(dnn_df)] = new_s

    dnn_df.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    gene_symbol_to_id()
