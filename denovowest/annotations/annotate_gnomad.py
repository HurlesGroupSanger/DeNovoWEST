#!/usr/bin/env python

import pandas as pd
import click
from pysam import VariantFile
from itertools import groupby, count
import numpy as np


def load_gnomad(gnomad, chrom, start, end):
    """
    Access to specific region in gnomad VCF file and builds a data frame
    Args:
        gnomad (pysam.TabixFile): CADD score file
        chrom (str): chromosome identifier
        start (int): beginning position of the region
        end (int): terminating position of the region
    """
    gnomad_data = []
    for record in gnomad.fetch(chrom, start, end):
        for i, alt in enumerate(record.alts):
            if "AF" in record.info and record.info["AF"][i] is not None:
                dic = {
                    "chrom": str(record.chrom),
                    "pos": int(record.pos),
                    "ref": str(record.ref),
                    "alt": str(alt),
                    "maf": float(record.info["AF"][i]),
                }
                gnomad_data.append(dic)
    return pd.DataFrame(gnomad_data)


def as_range(region):
    """
    Returns genomic region boundaries

    Args:
        region (list): genomic region positions

    Returns:
        tuple: beginning and terminating position of the genomic region
    """
    l = list(region)
    return l[0], l[-1]


@click.command()
@click.argument("rates")
@click.argument("gnomad")
@click.argument("output")
def annotate_gnomad(rates, gnomad, output):
    """Adds gnomad MAF to rates file

    Args:
        rates (str): Path to rates file
        gnomad (str): Path to gnomad file
        output (str): Path to output file (merged dataframe)
    """

    # Load rates file
    rates_df = pd.read_table(rates)

    # Load gnomad VCF
    gnomad_df = VariantFile(gnomad)

    list_merged_df = list()
    for gene_id, gene_rates_df in rates_df.groupby("symbol"):
        chrom = gene_rates_df.chrom.values[0]

        # Split each gene in contiguous block (i.e. exons) and load gnomad data
        # (memory and performance issue)
        list_block_df = list()
        for _, block in groupby(sorted(set(gene_rates_df["pos"])), key=lambda n, c=count(): n - next(c)):
            start, end = as_range(block)

            block_gnomad_df = load_gnomad(gnomad_df, chrom, start - 1, end)
            list_block_df.append(block_gnomad_df)

        gene_gnomad_df = pd.concat(list_block_df)

        # Merge rates with Gnomad data
        if gene_gnomad_df.empty:
            merged_gene_df = gene_rates_df
            gene_rates_df["maf"] = np.nan
        else:
            merged_gene_df = gene_rates_df.merge(gene_gnomad_df, how="left", on=["chrom", "pos", "ref", "alt"])
        list_merged_df.append(merged_gene_df)

    # Combine results for all genes
    merged_df = pd.concat(list_merged_df)
    merged_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    annotate_gnomad()
