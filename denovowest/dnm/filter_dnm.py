#!/usr/bin/env python
import click
import pandas as pd


@click.command()
@click.argument("dnm")
@click.argument("gene_list")
@click.option("--output")
def filter_dnm(dnm, gene_list, output):
    """_summary_

    Args:
        dnm (_type_): _description_
        gene_list (_type_): _description_
        output (_type_): _description_
    """
    dnm_df = load_dnm(dnm)
    genes = load_gene_list(gene_list)

    dnm_df = filter_dnm_df(dnm_df, genes)

    dnm_df.to_csv(output, sep="\t", index=False)


def filter_dnm_df(dnm_df, gene_list):
    """_summary_

    Args:
        dnm (_type_): _description_
        gene_list (_type_): _description_
    """

    dnm_df = dnm_df.loc[dnm_df["gene_id"].isin(gene_list)]

    return dnm_df


def load_dnm(dnm):
    """_summary_

    Args:
        dnm (_type_): _description_
    """

    dnm_df = pd.read_csv(dnm, sep="\t")

    dnm_df["gene_id"] = format_gene_id(list(dnm_df["gene_id"]))

    return dnm_df


def load_gene_list(gene_list):
    """_summary_

    Args:
        gene_list (_type_): _description_
    """

    with open(gene_list, "r") as f:
        genes = f.readlines()

    genes = format_gene_id(genes)

    return genes


def format_gene_id(gene_list):
    """_summary_

    Args:
        gene_list (_type_): _description_
    """
    gene_list = [gene_id.strip() for gene_id in gene_list]
    # gene_list = [gene_id.replace("gene:", "") for gene_id in gene_list]

    return gene_list


if __name__ == "__main__":
    test = filter_dnm()
