#!/usr/bin/env python
import pandas as pd
import numpy as np
import click
import json
import pysam
import subprocess
import utils
import logging


def rates_stats(rates_df):
    """
    Compute statistics such as missing annotations for each gene

    Args:
        rates_df (pd.DataFrame): rates
    """

    rates_stats_dict = dict()

    for gene_id, generates_df in rates_df.groupby("gene_id"):
        rates_stats_dict[gene_id] = gene_stats(gene_id, generates_df)

    return rates_stats_dict


def gene_stats(gene_id, generates_df):
    """
    Compute statistics such as missing annotations for each molecular consequence
    and each annotation.

    Args:
        gene_id (str) : gene identifier
        generates_df (pd.DataFrame): rates for a given gene
    """

    gene_dict = dict()

    for consequence, consequence_df in generates_df.groupby("consequence"):
        gene_dict[consequence] = consequence_stats(consequence, consequence_df)

    return gene_dict


def consequence_stats(consequence, consequence_df):
    """
    Compute statistics such as missing annotations for each molecular consequence
    and each annotation.

    Args:
        consequence (str): molecular consequence
        consequence_df (pd.DataFrame)
    """

    consequence_dict = dict()

    nb_mutations = consequence_df.shape[0]
    nb_missing_prob = consequence_df["prob"].isna().sum()
    consequence_dict["nb_mutations"] = int(nb_mutations)
    consequence_dict["nb_missing_prob"] = int(nb_missing_prob)

    scores = list(consequence_df.select_dtypes(include=["float"]).columns)

    for score in scores:
        score_dict = dict()
        nb_missing_score = consequence_df[score].isna().sum()
        deciles = consequence_df[score].quantile([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

        score_dict["nb_missing_score"] = int(nb_missing_score)
        if nb_missing_score != nb_mutations:
            score_dict["deciles"] = list(deciles)

        consequence_dict[score] = score_dict

    return consequence_dict


def handle_nas(rates_df):
    """
    Turns "." and empty cells into NaNs in order to convert scores columns to float

    Args:
        rates_df (pd.DataFrame): rates
    """

    # TODO : Better handle this, the rates file should probably has missing values rather than "."
    rates_df.replace({".": np.nan, "": np.nan}, inplace=True)
    for column in rates_df.columns:
        if column in ["pos", "prob"]:
            continue
        try:
            rates_df[column] = rates_df[column].astype(float)
        except ValueError:
            continue

    return rates_df


def warning_genes(rates_stats_dict):
    """
    Find genes with too many missing scores for each annotation
    and log them.

    Args:
        rates_stats_dict (dict): annotation statistics
    """

    THRESHOLD_PERCENTAGE_MISSING_ANNOTATIONS = 0.1

    logger = logging.getLogger("logger")
    for gene_id, gene_stats in rates_stats_dict.items():
        perc_missing_annotations = percentage_missing_annotations(gene_stats)

        for score, percentage in perc_missing_annotations.items():
            if percentage > THRESHOLD_PERCENTAGE_MISSING_ANNOTATIONS:
                logger.warning(f"{gene_id} | {score} | {round(percentage * 100,2)}% missing annotations")


def percentage_missing_annotations(gene_stats):
    """
    Compute the percentage of missing scores for non-synonymous mutations.

    Args:
        gene_stats (dict): annotation statistics at the gene level

    Returns:
        dict: percentage of missing scores per annotation
    """

    consequences_to_consider = ["missense", "stop_gained", "splice_acceptor", "splice_donor", "start_lost", "stop_lost"]

    nb_mutation = 0
    nb_missing_annotations = dict()
    for consequence, consequence_stats in gene_stats.items():
        if consequence in consequences_to_consider:
            nb_mutation += consequence_stats["nb_mutations"]
            for score, score_stats in consequence_stats.items():
                if not (isinstance(score_stats, dict)):
                    continue
                if score in nb_missing_annotations.keys():
                    nb_missing_annotations[score] += score_stats["nb_missing_score"]
                else:
                    nb_missing_annotations[score] = score_stats["nb_missing_score"]

    for score, count in nb_missing_annotations.items():
        nb_missing_annotations[score] = count / nb_mutation

    return nb_missing_annotations


@click.command()
@click.argument("rates")
@click.argument("output")
@click.option("--fullgenome", is_flag=True)
def main(rates, output, fullgenome):
    """
    This script provides informations on an annotated mutation rates file
    such as the percentage of functional variants missing a score for each
    annotation.

    It is mainly dedicated to be run in the nextflow pipeline on subsetted
    rates file but can also be used on a full genome indexed rates file by using
    the --fullgenome flag.

    Args:
        rates (str): mutation rates file
        output (str): JSON file storing per-gene statistics
        fullgenome (boolean) : to use with full genome indexed mutation rates file
    """

    utils.init_log()
    logger = logging.getLogger("logger")

    # If using a full genome indexed rates file
    if fullgenome:
        rates_file = pysam.TabixFile(rates)

        # Retrieve rates file header
        cmd = f"zcat {rates} | head -n 1"
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        rates_file_columns = res.stdout.strip().split("\t")

        # Load one chromosome at a time (for memory sake) and compute stats
        chromosomes = ["chr" + str(x) for x in list(range(1, 23))] + ["chrX"]
        rates_stats_dict = dict()
        for chromosome in chromosomes:
            try:
                rates_chrom_df = pd.DataFrame(
                    [line.split("\t") for line in rates_file.fetch(chromosome)], columns=rates_file_columns
                )
            except ValueError:
                logger.warning(f"# {chromosome} not found in index")
                continue

            rates_chrom_df = handle_nas(rates_chrom_df)
            rates_chrom_stats_dict = rates_stats(rates_chrom_df)

            rates_stats_dict.update(rates_chrom_stats_dict)

    # Otherwise if using a subsetted rates file (e.g. from nextflow pipeline)
    else:
        rates_df = pd.read_csv(
            rates,
            sep="\t",
            dtype={
                "gene_id": str,
                "chrom": str,
                "pos": int,
                "ref": str,
                "alt": str,
                "consequence": str,
                "prob": float,
            },
        )

        rates_df = handle_nas(rates_df)
        rates_stats_dict = rates_stats(rates_df)

    # Export stats in a JSON file
    with open(output, "w") as f:
        json.dump(rates_stats_dict, f, indent=4)

    # Log genes with too many missing annotations
    warning_genes(rates_stats_dict)


if __name__ == "__main__":
    main()
