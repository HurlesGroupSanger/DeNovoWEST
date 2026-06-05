#!/usr/bin/env python
import logging
import sys
from itertools import count, groupby

import click
import pandas as pd
import pysam

from denovowest.utils.io_helpers import as_range
from denovowest.utils.log import init_log


def load_cadd(cadd_file: pysam.TabixFile, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Fetch CADD scores for a genomic region.

    Args:
        cadd_file: CADD score tabix file.
        chrom: Chromosome identifier (without "chr" prefix).
        start: 0-based start position.
        end: End position.

    Returns:
        DataFrame with columns [chrom, pos, ref, alt, raw, score].
    """

    def _parse(line: str) -> dict:
        chrom, pos, ref, alt, raw, scaled = line.split("\t")
        return {
            "chrom": chrom,
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
            "raw": float(raw),
            "score": float(scaled),
        }

    return pd.DataFrame([_parse(x) for x in cadd_file.fetch(chrom, start, end)])


@click.command()
@click.argument("rates")
@click.argument("cadd")
@click.argument("output")
def annotate_cadd(rates: str, cadd: str, output: str) -> None:
    """Add CADD scores to a rates file.

    Args:
        rates: Path to the rates file.
        cadd: Path to the CADD tabix file.
        output: Path to the output file.
    """

    init_log()
    logger = logging.getLogger("logger")

    # Load rates file
    rates_df = pd.read_table(rates, dtype={"chrom": str, "pos": int, "ref": str, "alt": str})
    if rates_df.empty:
        logger.warning("Rates file is empty")
        rates_df["raw"] = pd.NA
        rates_df["score"] = pd.NA
        rates_df.to_csv(output, sep="\t", index=False)
        sys.exit(0)

    # Depending on the gff, chromosome can be defined as "chrX" or just "X"
    add_chr = str(rates_df.iloc[0].chrom).startswith("chr")

    # Load CADD tabix file
    cadd_tabix = pysam.TabixFile(cadd)

    # For each gene
    list_merged_df = []
    for gene_id, gene_rates_df in rates_df.groupby("gene_id"):
        chrom = str(gene_rates_df.chrom.values[0]).replace("chr", "")

        # Split each gene in contiguous blocks (i.e. exons) and load CADD scores
        list_block_df = []
        for _, block in groupby(sorted(set(gene_rates_df["pos"])), key=lambda n, c=count(): n - next(c)):
            start, end = as_range(block)
            try:
                block_cadd_df = load_cadd(cadd_tabix, chrom, start - 1, end)
                list_block_df.append(block_cadd_df)
            except ValueError:
                logger.warning(f"No CADD scores for {gene_id} at {chrom}:{start}-{end}")
                continue

        # Merge rates with CADD scores
        if list_block_df:
            gene_cadd_df = pd.concat(list_block_df)
            if add_chr:
                gene_cadd_df["chrom"] = "chr" + gene_cadd_df["chrom"]
            merged_gene_df = gene_rates_df.merge(gene_cadd_df, how="left", on=["chrom", "pos", "ref", "alt"])
        else:
            merged_gene_df = gene_rates_df.copy()
            merged_gene_df["raw"] = pd.NA
            merged_gene_df["score"] = pd.NA

        list_merged_df.append(merged_gene_df)

    # Combine results for all genes
    merged_df = pd.concat(list_merged_df)
    merged_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    annotate_cadd()
