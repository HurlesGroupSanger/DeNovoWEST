#!/usr/bin/env python
import click
import pandas as pd
import logging
import re
import gffutils

from denovowest.utils.log import init_log
from denovowest.utils.params import CDS_OFFSET
from denovowest.utils.io_helpers import load_gff


@click.command()
@click.argument("dnm")
@click.argument("gene_list")
@click.option("--output_kept_dnm")
@click.option("--output_discarded_dnm")
@click.option("--gff")
def filter_dnm(dnm, gene_list, output_kept_dnm, output_discarded_dnm, gff):
    """
    Remove from DNM table the variants in genes not found in gene list.
    The gene list has been either provided by the user or built from the rates file.

    If the GFF option is provided, variants not found in CDS regions will be discarded.

    Args:
        dnm (str): path to DNM file
        gene_list (str): path to gene list file
        output_kept_dnm (str): path to kept DNM table
        output_discarded_dnm (str): path to filtered DNM table
        gff (str): path to GFF file
    """

    genes = load_gene_list(gene_list)
    dnm_df = load_dnm(dnm, genes)

    # Filter DNM table on gene list
    dnm_df, dnm_discarded_df = filter_on_gene_list(dnm_df, genes)

    # Filter DNM table on CDS regions
    if gff:
        dnm_df, dnm_discarded_gff_df = filter_on_gff(dnm_df, gff)
        dnm_discarded_df = pd.concat([dnm_discarded_df, dnm_discarded_gff_df])

    # Export DNM tables
    dnm_df.to_csv(output_kept_dnm, sep="\t", index=False)
    dnm_discarded_df.to_csv(output_discarded_dnm, sep="\t", index=False)

    # Log stats
    log_stats(dnm_df, dnm_discarded_df, output_discarded_dnm)


def filter_on_gene_list(dnm_df, gene_list):
    """
    Remove from DNM table the variants in genes not found in gene list

    Args:
        dnm_df (pd.DataFrame): DNM table
        gene_list (list): gene list
    """

    dnm_kept_df = dnm_df.loc[dnm_df["gene_id"].isin(gene_list)]
    dnm_discarded_df = dnm_df.loc[~dnm_df["gene_id"].isin(gene_list)].copy()
    dnm_discarded_df["reason"] = "not_in_gene_list"

    return dnm_kept_df, dnm_discarded_df


def filter_on_gff(dnm_df, gff):
    """
    Remove from DNM table the variants not found in CDS regions

    Args:
        dnm_df (pd.DataFrame): DNM table
        gff (str): path to GFF file or gffutils database
    """

    logger = logging.getLogger("logger")

    gff_db = load_gff(gff)

    idx_discarded = list()
    reason_discarded = list()
    for idx, dnm in dnm_df.iterrows():

        pos = dnm["pos"]
        gene_id = dnm["gene_id"]

        try:
            gene = gff_db[gene_id]
        except gffutils.exceptions.FeatureNotFoundError:

            # Rare case where a user provide a GFF that differs
            # from the gene used to annotate the DNM
            logger.warning(f"{gene_id} not in GFF")
            idx_discarded.append(idx)
            reason_discarded.append("gene_not_in_gff")
            continue

        if not is_in_cds(gff_db, gene, pos):
            idx_discarded.append(idx)
            reason_discarded.append("gene_not_in_cds")

    dnm_discarded_df = dnm_df.loc[idx_discarded].copy()
    dnm_discarded_df["reason"] = reason_discarded
    dnm_kept_df = dnm_df.drop(idx_discarded)

    return dnm_kept_df, dnm_discarded_df


def is_in_cds(gff_db, gene, pos):
    """
    Check if a given DNM is in a CDS region defined in the GFF file

    Args:
        gff_db (gffutils.FeatureDB): gffutils database
        gene (gffutils.Feature): gene in GFFDB
        pos (int) : genomic loci of the variant
    """

    for transcript in gff_db.children(gene, level=1):
        for cds in gff_db.children(transcript, featuretype="CDS", order_by="start"):
            if (pos >= cds.start - CDS_OFFSET) and (pos <= cds.end + CDS_OFFSET):
                return True

    return False


def load_dnm(dnm, genes):
    """
    Load DNM file

    Args:
        dnm (str): path to DNM file
        genes (list) : list of gene identifiers
    """

    dnm_df = pd.read_csv(dnm, sep="\t")
    dnm_df["gene_id"] = format_gene_id(list(dnm_df["gene_id"]), genes)

    return dnm_df


def load_gene_list(gene_list):
    """
    Load gene list that will be used to filter the DNM table

    Args:
        gene_list (str): path to gene list
    """

    with open(gene_list, "r") as f:
        genes = f.readlines()

    genes = [gene_id.strip() for gene_id in genes]

    return genes


def format_gene_id(genes_in_dnm, genes_list):
    """
    Get rid of trailing spaces after gene ids and matches gene identifiers from rates file

    Args:
        genes_in_dnm (list): list of genes in DNM file
        genes_list (list) : list of genes provided by the user or taken from the rates file
    """

    # Strip genes identifiers to remove typo mistakes
    genes_in_dnm = [gene_id.strip() for gene_id in genes_in_dnm]

    # If the DNM file gene identifiers do not contain version number but the rates file does, we suffix the gene identifier with the version
    pattern_ensembl_version = r"^ENSG\d{11}\.\d+$"
    pattern_ensembl_id_only = r"^ENSG\d{11}$"
    if bool(re.match(pattern_ensembl_id_only, genes_in_dnm[0])) & bool(
        re.match(pattern_ensembl_version, genes_list[0])
    ):
        genes_dict = {x.split(".")[0]: x for x in genes_list}
        genes_in_dnm = [genes_dict[x] if x in genes_dict.keys() else x for x in genes_in_dnm]

    return genes_in_dnm


def log_stats(dnm_kept_df, dnm_discarded_df, output_discarded_dnm):
    """
    Log stats about the DNM table filtering

    Args:
        dnm_kept_df (pd.DataFrame): DNM kept
        dnm_discarded_df (pd.DataFrame): DNM filtered out
        output_discarded_dnm (str): path to filtered DNM table
    """
    nb_dnm = dnm_kept_df.shape[0] + dnm_discarded_df.shape[0]
    logger = logging.getLogger("logger")

    if dnm_discarded_df.empty:
        logger.info(f"All DNM ({nb_dnm}) were retained")
    else:
        nb_dnm_discarded = dnm_discarded_df.shape[0]
        logger.warning(f"{nb_dnm_discarded}/{nb_dnm} DNM were discarded, of which :")
        for reason, count in dict(dnm_discarded_df.reason.value_counts()).items():
            logger.warning(f"- {reason} : {count} DNM")
        logger.warning(f"Check filtered DNM table : {output_discarded_dnm}")


if __name__ == "__main__":
    init_log()
    filter_dnm()
