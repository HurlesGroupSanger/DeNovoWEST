#!/usr/bin/env python
import pandas as pd
import click
import pysam
import sys
from itertools import groupby, count

is_chr_prefixed_variants = False
is_chr_prefixed_annotation = False


def is_chr_prefixed(df_path):

    df = pd.read_csv(df_path, sep="\t", nrows=100)

    for key in ["chrom", "chr", "#chr", "#chrom"]:
        if key in df.columns:
            is_chr_prefixed = str(df[key].iloc[0]).startswith("chr")
            break

    # TODO : Better handle this case
    if "is_chr_prefixed" not in locals():
        is_chr_prefixed = str(df.iloc[0, 0]).startswith("chr")

    return is_chr_prefixed


def load_variants_file(variants):
    """_summary_

    Args:
        variants (_type_): _description_
    """

    variants_df = pd.read_table(variants, dtype={"chrom": str, "pos": int, "ref": str, "alt": str})

    return variants_df


def annotate(variants_df, annotation_df, columns):

    columns_indices, columns_names = extract_columns_indices(annotation_df, columns)

    list_annotated_df = list()
    for gene_id, gene_df in variants_df.groupby("gene_id"):

        annotated_gene_df = annotate_gene(gene_id, gene_df, annotation_df, columns_indices)
        list_annotated_df.append(annotated_gene_df)

    # Combine results for all genes
    annotated_df = pd.concat(list_annotated_df)
    annotated_df.rename(dict(zip([str(x) for x in columns_indices], columns_names)), axis=1, inplace=True)

    return annotated_df


def extract_columns_indices(annotation_df, columns):
    """
    Extract the indices of the column to use from the annotation file

    Args:
        annotation_df (pysam.TabixFile): annotation indexed file
        columns (str): columns to extract from the annotation file
    """

    all_columns = annotation_df.header[0].split("\t")

    # If no columns have been specified we annotate with all non standard columns
    if columns:
        columns_to_extract = columns.split(",")
    else:
        columns_to_extract = all_columns[4:]

    # Retrieve indices of columns to extract
    columns_indices = list()
    columns_names = list()
    for idx, name in enumerate(all_columns):
        if name in columns_to_extract:
            columns_indices.append(idx)
            columns_names.append(name)

    # If a score can not be retrieved from dbNSFP we stop here
    columns_not_found = set(columns_to_extract) - set(columns_names)
    if columns_not_found:
        print(f"ERROR : The following columns could not be retrieved from the annotation file : {columns_not_found}")
        sys.exit(1)

    return columns_indices, columns_names


def annotate_gene(gene_id, gene_df, annotation_df, columns):

    gene_annotation_df = retrieve_annotation(gene_df, annotation_df, columns)

    if not gene_annotation_df.empty:
        gene_annotated_df = gene_df.merge(gene_annotation_df, how="left", on=["chrom", "pos", "ref", "alt"])
    else:
        print(f"WARNING : No annotations for gene {gene_id}")
        gene_annotated_df = gene_df

    # for column in columns:
    #     if column not in gene_annotated_df.columns:
    #         gene_annotated_df[column] = "."

    return gene_annotated_df


def retrieve_annotation(gene_df, annotation_df, columns):
    """_summary_

    Args:
        gene_df (_type_): _description_
        annotation_df (_type_): _description_

    Returns:
        _type_: _description_
    """

    gene_chrom = str(gene_df.chrom.values[0]).replace("chr", "")
    # Split each gene in contiguous block (i.e. exons) and retrieve annotations
    list_block_df = list()
    for _, block in groupby(sorted(set(gene_df["pos"])), key=lambda n, c=count(): n - next(c)):

        # Get the coordinates of the current block
        start, end = as_range(block)

        # Retrieve annotation for the current block
        try:
            block_cadd_df = load_annotation(annotation_df, gene_chrom, start - 1, end, columns)
            list_block_df.append(block_cadd_df)
        except ValueError:
            continue

    # Concatenate into gene wide annotation
    if list_block_df:
        gene_annotations_df = pd.concat(list_block_df)

        # Add the chr prefix if found in the variants file
        if not gene_annotations_df.empty and is_chr_prefixed_variants:
            gene_annotations_df["chrom"] = ["chr" + x for x in gene_annotations_df["chrom"]]
    else:
        gene_annotations_df = pd.DataFrame()

    return gene_annotations_df


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


def load_annotation(annotation_file, chrom, start, end, columns):
    """
    Access to specific region in annotation file
    Args:
        annotation_file (pysam.TabixFile): annotation file
        chrom (str): chromosome identifier
        start (int): beginning position of the region
        end (int): terminating position of the region
    """

    list_records = list()
    for record in annotation_file.fetch(chrom, start, end, columns):
        list_records.append(parse(record, columns))

    records_df = pd.DataFrame(list_records)

    # TODO : handle annotation duplication (especially in the context of gene/transcript id)
    records_df.drop_duplicates(inplace=True)

    return records_df


def parse(line, columns):

    record = line.split("\t")

    record_dict = dict()
    record_dict["chrom"] = record[0]
    record_dict["pos"] = int(record[1])
    record_dict["ref"] = record[2]
    record_dict["alt"] = record[3]

    for column in columns:
        record_dict[str(column)] = record[column]

    return record_dict


@click.command()
@click.argument("variants", type=click.Path(exists=True))
@click.argument("annotation", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
@click.option("--columns", type=str, help="Columns to use in the annotation file. Should be comma separated")
def cli(variants, annotation, output, columns):
    """
    Annotate a TSV variants file (DNM or rates) with information coming from a
    TAB-delimited genome position file indexed with tabix
    (https://www.htslib.org/doc/tabix.html)

    Args:
        variants (str): variants (DNM, rates) file
        annotation (str): annotation file
        output (str): variants file annotated
        columns (str): columns to use in the annotation file
    """

    global is_chr_prefixed_variants, is_chr_prefixed_annotation
    is_chr_prefixed_variants = is_chr_prefixed(variants)
    is_chr_prefixed_annotation = is_chr_prefixed(annotation)

    # Load variants file
    variants_df = load_variants_file(variants)

    # Load indexed annotation file
    annotation_df = pysam.TabixFile(annotation)

    # Annotate variants file
    annotated_df = annotate(variants_df, annotation_df, columns)

    # Export annotated file
    annotated_df.to_csv(output, sep="\t", index=False, na_rep=".")


if __name__ == "__main__":
    cli()
