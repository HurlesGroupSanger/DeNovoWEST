#!/usr/bin/env python
import click
import pandas as pd
import logging
import re
import gffutils
import numpy as np

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

    # Build gene CDS intervals
    gene_ids = list(dnm_df["gene_id"].unique())
    genes_cds = build_gene_cds_intervals(gff_db, gene_ids, cds_offset=CDS_OFFSET)
    logger.info("Built gene CDS intervals from GFF")

    # Loop over genes
    list_new_dnm_df = []
    for gene_id, gene_dnm_df in dnm_df.groupby("gene_id"):

        # If gene not in GFF, all its DNM are filtered out
        gene_cds = genes_cds[gene_id]
        if gene_cds is None:
            logger.warning(f"{gene_id} not in GFF")
            gene_dnm_df["in_cds"] = False
            gene_dnm_df["reason"] = "gene_not_in_gff"
            list_new_dnm_df.append(gene_dnm_df)
            continue

        # Loop over DNM in the gene to check if they are in CDS regions
        hit = np.zeros(gene_dnm_df.shape[0], dtype=bool)
        for cds_start, cds_end in gene_cds:
            hit |= (gene_dnm_df.pos >= cds_start) & (gene_dnm_df.pos <= cds_end)

        gene_dnm_df["in_cds"] = hit
        gene_dnm_df["reason"] = gene_dnm_df["in_cds"].apply(lambda x: "not_in_cds" if not x else "")
        list_new_dnm_df.append(gene_dnm_df)

    # Combine all DNM
    new_dnm_df = pd.concat(list_new_dnm_df)
    kept_df = new_dnm_df.loc[new_dnm_df["in_cds"]].copy()
    discarded_df = new_dnm_df.loc[~new_dnm_df["in_cds"]].copy()

    return kept_df.drop(["in_cds", "reason"], axis=1), discarded_df.drop(["in_cds"], axis=1)


def build_gene_cds_intervals(gff_db, gene_ids, cds_offset):
    """
    Build merged CDS intervals for each gene in gene_ids.

    Returns mapping: gene_id -> list of (start, end) sorted, non-overlapping.
    If a gene is missing in the GFF, maps to None.
    """
    intervals = {}

    for gid in gene_ids:
        try:
            gene = gff_db[gid]
        except gffutils.exceptions.FeatureNotFoundError:
            intervals[gid] = None
            continue

        raw_intervals = []
        for transcript in gff_db.children(gene, level=1):
            for cds in gff_db.children(transcript, featuretype="CDS", order_by="start"):
                start = max(1, cds.start - cds_offset)
                end = cds.end + cds_offset
                raw_intervals.append((start, end))

        if not raw_intervals:
            intervals[gid] = []
            continue

        intervals[gid] = raw_intervals

    return intervals


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
