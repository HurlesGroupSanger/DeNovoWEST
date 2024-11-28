#!/usr/bin/env python
import pandas as pd
import json
import click
import os


@click.command()
@click.argument("json_stats")
@click.argument("outdir")
def main(json_stats, outdir):
    """
    This script loads a per-gene statistics JSON file containing informations
    (such as the number of missing annotation for a given CEP on a given gene
    and for a given type of molecular consequence) and creates data frame that summarize
    these statistics.

    Args:
        json_stats (str): per-gene statistics JSON file
        outdir (str): output directory
    """

    # Load JSON file containing statistics per gene
    with open(json_stats) as f:
        stats = json.load(f)

    # Create the output directory if not already existing
    os.makedirs(outdir, exist_ok=True)

    # Retrieve the list of all molecular consequences inferred by bcftools csq and the CEP used
    # to assign scores to every mutation in the rates file
    consequences = list_consequences(stats)
    ceps = list_ceps(stats)

    # Calculate the coverage for each gene for each CEP for each molecular consequence
    for consequence in consequences:
        coverage_per_consequence_df = coverage_per_consequence(stats, consequence, ceps)
        coverage_per_consequence_df.to_csv(f"{outdir}/coverage_{consequence.replace('&', '-')}.tsv", sep="\t")


def list_consequences(stats):
    """
    Return the list of bcftoolscsq consequences assigned to every locus
    in the rates file.

    Args:
        stats (dict): Per-gene statistics
    """

    list_consequences = list()
    for gene_id, gene_dict in stats.items():
        list_consequences += list(gene_dict.keys())

    list_consequences = list(set(list_consequences))

    return list_consequences


def list_ceps(stats):
    """
    Returns the list of Computational Effect Predictors (CEP) used to
    assign a score to every locus in the rates file

    Args:
        stats (dict): Per-gene statistics
    """

    list_ceps = list()
    for gene_id, gene_dict in stats.items():
        list_ceps = list(gene_dict["missense"].keys())
        break

    list_ceps = list(set(list_ceps) - {"nb_mutations", "nb_missing_prob", "prob", "chrom"})

    return list_ceps


def coverage_per_consequence(stats, consequence, ceps):
    """
    For each CEP, retrieves the percentage of loci annotated in a given
    gene for a given consequence

    Args:
        stats (dict): Per-gene statistics
        consequence (str): Molecular consequence (e.g. missense)
        ceps (list): List of CEP (e.g. AlphaMissense_score)
    """

    consequence_dict = {
        gene_id: gene_dict[consequence] if consequence in gene_dict.keys() else {}
        for (gene_id, gene_dict) in stats.items()
    }

    return compute_percentage_annotated(consequence_dict, ceps)


def compute_percentage_annotated(consequence_dict, ceps):
    """
    For each CEP, retrieves the percentage of loci annotated in a given
    gene for a given consequence

    Args:
        consequence_dict (dict): statistics for a given gene and a given type of molecular consequence
        ceps (list): list of CEP used to score the variants
    """

    percentage_annotated = dict()
    for gene_id, gene_data in consequence_dict.items():

        # If the current gene do not have any variants of the current consequence type we skip it
        if not gene_data:
            continue

        # Otherwise we calculate the percentage of loci annotated with every CEP used to annotate the rates file
        nb_mutations = gene_data["nb_mutations"]
        percentage_annotated[gene_id] = dict()
        for cep in ceps:
            percentage_annotated[gene_id][cep] = round(gene_data[cep]["nb_missing_score"] * 100 / nb_mutations, 2)

    return pd.DataFrame(percentage_annotated).transpose()


if __name__ == "__main__":
    main()
