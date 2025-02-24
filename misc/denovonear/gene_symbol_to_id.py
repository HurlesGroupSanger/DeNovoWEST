#!/usr/bin/env python
import click
import pandas as pd
import gffutils


@click.command()
@click.argument("dnn_results", type=click.Path(exists=True))
@click.argument("gff", type=click.Path(exists=True))
@click.argument("out", type=str)
def gene_symbol_to_id(dnn_results, gff, out):
    """_summary_

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
    for gene in gff_db.features_of_type("gene", order_by="start"):
        gene_dict[gene["Name"][0]] = gene["ID"][0]

    # Get the gene identifier
    dnn_df["gene_id"] = dnn_df["gene_symbol"].map(gene_dict)

    # Export Denovonear with ensembl IDs
    dnn_df = dnn_df[["gene_id", "gene_symbol", "mutation_category", "events_n", "dist", "probability"]]
    dnn_df.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    gene_symbol_to_id()
