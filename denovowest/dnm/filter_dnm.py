#!/usr/bin/env python
import click
import pandas as pd


@click.command()
@click.argument("dnm")
@click.argument("gene_list")
@click.option("--output_kept_dnm")
@click.option("--output_discarded_dnm")
def filter_dnm(dnm, gene_list, output_kept_dnm, output_discarded_dnm):
    """Remove from DNM table the variants in genes not found in gene list.
    The gene list has been either provided by the user or built from the rates file.

    Args:
        dnm (str): path to DNM file
        gene_list (str): path to gene list file
        output_kept_dnm (str): path to kept DNM table
        output_discarded_dnm (str): path to filtered DNM table
    """

    dnm_df = load_dnm(dnm)
    nb_dnm = dnm_df.shape[0]
    genes = load_gene_list(gene_list)

    dnm_df, dnm_discarded_df = filter_dnm_df(dnm_df, genes)

    dnm_df.to_csv(output_kept_dnm, sep="\t", index=False)
    dnm_discarded_df.to_csv(output_discarded_dnm, sep="\t", index=False)

    if dnm_discarded_df.empty:
        print(f"All DNM ({nb_dnm}) were retained.")
    else:
        print(
            f"{dnm_discarded_df.shape[0]}/{nb_dnm} DNM in {dnm_discarded_df['gene_id'].nunique()} genes have been discarded. See {output_discarded_dnm} for more details."
        )


def filter_dnm_df(dnm_df, gene_list):
    """Remove from DNM table the variants in genes not found in gene list

    Args:
        dnm_df (pd.DataFrame): DNM table
        gene_list (list): gene list
    """

    dnm_kept_df = dnm_df.loc[dnm_df["gene_id"].isin(gene_list)]
    dnm_discarded_df = dnm_df.loc[~dnm_df["gene_id"].isin(gene_list)]

    return dnm_kept_df, dnm_discarded_df


def load_dnm(dnm):
    """Load DNM file

    Args:
        dnm (str): path to DNM file
    """

    dnm_df = pd.read_csv(dnm, sep="\t")
    dnm_df["gene_id"] = format_gene_id(list(dnm_df["gene_id"]))

    return dnm_df


def load_gene_list(gene_list):
    """Load gene list that will be used to filter the DNM table

    Args:
        gene_list (str): path to gene list
    """

    with open(gene_list, "r") as f:
        genes = f.readlines()

    genes = format_gene_id(genes)

    return genes


def format_gene_id(gene_list):
    """
    Get rid of trailing spaces after gene ids

    Args:
        gene_list (list): list of genes in DNM file
    """
    gene_list = [gene_id.strip() for gene_id in gene_list]
    # TODO : Handle B37 GFF case where gene ids are prefixed by "gene:"
    # gene_list = [gene_id.replace("gene:", "") for gene_id in gene_list]

    return gene_list


if __name__ == "__main__":
    test = filter_dnm()
