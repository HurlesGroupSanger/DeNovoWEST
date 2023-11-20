#!/usr/bin/env python
import pandas as pd
import click
import pysam
import sys
from itertools import groupby, count


def load_cadd(cadd_file, chrom, start, end):
    """
    Access to specific region in CADD file.
    Args:
        cadd_file (pysam.TabixFile): CADD score file
        chrom (str): chromosome identifier
        start (int): beginning position of the region
        end (int): terminating position of the region
    """

    def parse(line):
        chrom, pos, ref, alt, raw, scaled = line.split("\t")

        return {
            "chrom": chrom,
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
            "raw": float(raw),
            "score": float(scaled),
        }

    return pd.DataFrame([parse(x) for x in cadd_file.fetch(chrom, start, end)])


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
@click.argument("cadd")
@click.argument("output")
def annotate_cadd(rates, cadd, output):
    """Adds CADD score to rates file

    Args:
        rates (str): Path to rates file
        cadd (str): Path to CADD file
        output (str): Path to output file (merged dataframe)
    """

    # Load rates file
    rates_df = pd.read_table(rates, dtype={"chrom": str, "pos": int, "ref": str, "alt": str})
    if rates_df.empty:
        print("Rates file is empty")
        rates_df["raw"] = None
        rates_df["score"] = None
        rates_df.to_csv(output, sep="\t", index=False)
        sys.exit(0)

    # Depending on the gff, chromosome can be defined as "chrX" or just "X"
    if str(rates_df.iloc[0].chrom).startswith("chr"):
        add_chr = True
    else:
        add_chr = False

    # Load cadd file
    cadd_df = pysam.TabixFile(cadd)

    # For each gene
    list_merged_df = list()
    for gene_id, gene_rates_df in rates_df.groupby("gene_id"):
        chrom = str(gene_rates_df.chrom.values[0]).replace("chr", "")

        # Split each gene in contiguous block (i.e. exons) and load CADD scores
        # (memory and performance issue)
        list_block_df = list()
        for _, block in groupby(sorted(set(gene_rates_df["pos"])), key=lambda n, c=count(): n - next(c)):
            start, end = as_range(block)
            try:
                block_cadd_df = load_cadd(cadd_df, chrom, start - 1, end)
                list_block_df.append(block_cadd_df)
            except ValueError:
                continue

        # Merge rates with CADD score
        if list_block_df:
            gene_cadd_df = pd.concat(list_block_df)
            if add_chr:
                gene_cadd_df.chrom = "chr" + gene_cadd_df.chrom

            merged_gene_df = gene_rates_df.merge(gene_cadd_df, how="left", on=["chrom", "pos", "ref", "alt"])
        else:
            merged_gene_df = gene_rates_df

        list_merged_df.append(merged_gene_df)

    # Combine results for all genes
    merged_df = pd.concat(list_merged_df)
    merged_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    merged_df = annotate_cadd()
