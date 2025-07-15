#!/usr/bin/env python
import pandas as pd
import click
import sys
import logging

from denovowest.utils.log import init_log


@click.command()
@click.argument("rates")
@click.argument("shet")
@click.argument("output")
@click.option("--threshold", default=0.15, help="shet threshold to determine if gene is high or low")
def annotate_shet(rates, shet, output, threshold):
    """Adds shet boolean to rates file

    Args:
        rates (str): Path to rates file
        shet (str): Path to shet file
        output (str): Path to output file (merged dataframe)
    """

    init_log()
    logger = logging.getLogger("logger")

    # Load rates file
    rates_df = pd.read_table(rates)
    if rates_df.empty:
        logger.warning("Rates file is empty")
        rates_df["shethigh"] = None
        rates_df.to_csv(output, sep="\t", index=False)
        sys.exit(0)
    rates_df["shethigh"] = False

    # Load shet file
    shet_df = pd.read_table(shet, index_col=1)

    # For each gene
    nb_gene_not_in_shet = 0
    for gene_id_orig, gene_rates_df in rates_df.groupby("gene_id"):
        # Remove gene version if any
        gene_id = gene_id_orig.split(".")[0]
        gene_id = gene_id.replace("gene:", "")

        try:
            mean_shet = shet_df.loc[gene_id, "mean_s_het"]
            if mean_shet > threshold:
                rates_df.loc[rates_df["gene_id"] == gene_id_orig, "shethigh"] = True

        except KeyError:
            nb_gene_not_in_shet += 1

    logger.info(f"{nb_gene_not_in_shet}/{len(set(rates_df['gene_id']))} genes from rates file not found in shet")

    # Export results
    rates_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    merged_df = annotate_shet()
