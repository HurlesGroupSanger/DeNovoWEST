#!/usr/bin/env python
import pandas as pd
import click
import pysam
import sys
from itertools import groupby, count


def load_dbnsfp(dbnsfp_file, chrom, start, end, columns_indices):
    """
    Access to specific region in an dbNSFP indexed file
    Args:
        annotation_file (pysam.TabixFile): Annotation file containing CEP scores
        chrom (str): chromosome identifier
        start (int): beginning position of the region
        end (int): terminating position of the region
        columns_indices(list) : indices of the dbNSFP columns to retrieve
    """

    def parse(record):
        """
        Parse a dbNSFP record (line) and extract the data of interest

        Args:
            record (str): line in dbNSFP

        Returns:
            dict: informations about a mutation (pos, CEP scores)
        """

        record = record.split("\t")

        record_dict = dict()
        record_dict["chrom"] = record[0]
        record_dict["pos"] = int(record[1])
        record_dict["ref"] = record[2]
        record_dict["alt"] = record[3]

        for column_index in columns_indices:
            record_dict[str(column_index)] = record[column_index]

        return record_dict

    return pd.DataFrame([parse(x) for x in dbnsfp_file.fetch(chrom, start, end)])


def as_range(region):
    """
    Returns genomic region boundaries

    Args:
        region (list): genomic region positions

    Returns:
        tuple: beginning and terminating position of the genomic region
    """

    return region[0], region[-1]


def extract_columns_from_dbNFP(dbnsfp, annotation_names):
    """
    Extract the indices of the column to use in dbnsfp

    Args:
        dbnsfp (pysam.TabixFile): dbNSFP indexed file
        annotation_names (str): path to a file listing column names from dbNSFP to retrieve
    """

    with open(annotation_names, "r") as f:
        columns_to_extract = [x.strip() for x in f.readlines()]

    dbnsfp_columns = dbnsfp.header[0].split("\t")

    list_indices = list()
    found_column = list()
    for idx, name in enumerate(dbnsfp_columns):
        if name in columns_to_extract:
            list_indices.append(idx)
            found_column.append(name)

    # If a score can not be retrieved from dbNSFP we stop here
    columns_not_found = set(columns_to_extract) - set(found_column)
    if columns_not_found:
        print(f"ERROR : The following columns could not be retrieved from dbNFSP : {columns_not_found}")
        sys.exit(1)

    return list_indices, columns_to_extract


def load_rates_file(rates, output):
    """
    Load rates file

    Args:
        rates (_type_): _description_
        output (_type_): _description_

    Returns:
        _type_: _description_
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

    return rates_df, add_chr


@click.command()
@click.argument("rates")
@click.argument("dbnsfp")
@click.argument("annotation_names")
@click.argument("output")
def annotate_dbnsfp(rates, dbnsfp, annotation_names, output):
    """Adds CADD score to rates file

    Args:
        rates (str): Path to rates file
        dbnsfp (str): Path to dbNSFP file
        annotation_names (str): Path to a file listing annotations to extract and merge
        output (str): Path to output file (merged dataframe)
    """

    # Load rates file
    rates_df, add_chr = load_rates_file(rates, output)
    rates_df_columns = list(rates_df.columns)

    # Load dbnsfp file
    dbnsfp_df = pysam.TabixFile(dbnsfp)

    # Retrieve columns to extract from dbnsfp_df
    dbnsfp_columns_indices, dbnsfp_columns_names = extract_columns_from_dbNFP(dbnsfp_df, annotation_names)

    # For each gene
    list_merged_df = list()
    for gene_id, gene_rates_df in rates_df.groupby("gene_id"):
        chrom = str(gene_rates_df.chrom.values[0]).replace("chr", "")

        # Split each gene in contiguous block (i.e. exons) and load dbNSFP scores
        # (memory and performance issue)
        list_block_df = list()
        for _, block in groupby(sorted(set(gene_rates_df["pos"])), key=lambda n, c=count(): n - next(c)):
            start, end = as_range(list(block))
            try:
                block_dbnsfp = load_dbnsfp(dbnsfp_df, chrom, start - 1, end, dbnsfp_columns_indices)
                if not block_dbnsfp.empty:
                    list_block_df.append(block_dbnsfp)
            except ValueError:
                continue

        # Merge rates with CADD score
        if list_block_df:
            gene_dbnsfp_df = pd.concat(list_block_df)
            if add_chr:
                gene_dbnsfp_df.chrom = "chr" + gene_dbnsfp_df.chrom

            merged_gene_df = gene_rates_df.merge(gene_dbnsfp_df, how="left", on=["chrom", "pos", "ref", "alt"])
        else:
            merged_gene_df = gene_rates_df

        list_merged_df.append(merged_gene_df)

    # Combine results for all genes
    merged_df = pd.concat(list_merged_df)
    merged_df.columns = rates_df_columns + dbnsfp_columns_names
    merged_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    merged_df = annotate_dbnsfp()
