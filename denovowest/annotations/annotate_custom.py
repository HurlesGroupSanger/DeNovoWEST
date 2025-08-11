#!/usr/bin/env python
import pandas as pd
import click
import pysam
import sys
from itertools import groupby, count
import logging
import gzip

from denovowest.utils.log import init_log


is_chr_prefixed_variants = False
is_chr_prefixed_annotation = False

# TODO : Handle annotation with matching on gene identifier/symbol


def is_chr_prefixed(df_path):
    """
    Chromosomes can be defined as "N" or "chrN"
    If it differs between the input file and the annotation file it prevents records matching

    Args:
        df_path (str): variants or annotation file

    Returns:
        bool: whether the chromosomes are prefixed or not in {df_path}
    """

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
    """
    Load variants file (e.g. DNM/rates file)

    Args:
        variants (str): variants file
    """

    variants_df = pd.read_table(variants, dtype={"chrom": str, "pos": int, "ref": str, "alt": str})

    return variants_df


def annotate(variants_df, annotation_df, columns_indices, columns_names):
    """
    Annotate variants file with the selected columns from the annotation file.

    Args:
        variants_df (pd.DataFrame): variants
        annotation_df (pd.DataFrame): annotations
        columns_indices (list):  index of annotations to extract
        columns_names (list):  names of annotations to extract
    """

    list_annotated_df = list()
    for gene_id, gene_df in variants_df.groupby("gene_id"):

        annotated_gene_df = annotate_gene(gene_id, gene_df, annotation_df, columns_indices)
        list_annotated_df.append(annotated_gene_df)

    # Combine results for all genes
    annotated_df = pd.concat(list_annotated_df)
    annotated_df.rename(dict(zip([str(x) for x in columns_indices], columns_names)), axis=1, inplace=True)

    return annotated_df


def annotate_gene(gene_id, gene_df, annotation_df, columns):
    """
    Annotate the current gene

    Args:
        gene_id (str): current gene identifier
        gene_df (pd.DataFrame): variants found in the current gene
        annotation_df (pd.DataFrame): annotations
        columns (list):  indices of annotation to retrieve

    """

    logger = logging.getLogger("logger")

    gene_annotation_df = retrieve_annotation(gene_df, annotation_df, columns, gene_id)

    if not gene_annotation_df.empty:
        gene_annotated_df = gene_df.merge(gene_annotation_df, how="left", on=["chrom", "pos", "ref", "alt"])
    else:
        logger.warning(f"No annotations for gene {gene_id}")
        gene_annotated_df = gene_df

    assert gene_df.shape[0] == gene_annotated_df.shape[0]

    # for column in columns:
    #     if column not in gene_annotated_df.columns:
    #         gene_annotated_df[column] = "."

    return gene_annotated_df


def retrieve_annotation(gene_df, annotation_df, columns, gene_id):
    """
    Retrieve annoations for the current gene

    Args:
        gene_df (pd.DataFrame): _description_
        annotation_df (pd.DataFrame): _description_
        columns (list):  indices of annotation to retrieve
        gene_id (str) : gene identifier

    """

    logger = logging.getLogger("logger")

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

    # Remove duplicated rows if any
    if not gene_annotations_df.empty:
        nb_rows_before_duplication = gene_annotations_df.shape[0]
        gene_annotations_df = gene_annotations_df.loc[~gene_annotations_df[["chrom", "pos", "ref", "alt"]].duplicated()]
        if nb_rows_before_duplication != gene_annotations_df.shape[0]:
            logger.warning(f"There were multiple annotation available for {gene_id}")

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


def extract_columns_indices(annotation_df, columns, columns_file):
    """
    Extract the indices of the column to use from the annotation file

    Args:
        annotation_df (pysam.TabixFile): annotation indexed file
        columns (str): columns to extract from the annotation file, comma separated
        columns_file (str)
    """

    logger = logging.getLogger("logger")

    try:
        all_columns = annotation_df.header[0].split("\t")
    except IndexError:
        # If the annotation file header does not start with #, the Tabix header will be empty
        # and we need to retrieve it another way
        with gzip.open(annotation_df.filename, "rt") as f:
            all_columns = f.readline().strip().split("\t")

    # If the user provided the columns to extract in a separate file we read it
    if columns_file:
        columns_to_extract = utils.read_columns_from_file(columns_file)
    # If he provided them as a string we split it
    elif columns:
        columns_to_extract = columns.split(",")
    # If he did not provide any column we retrieve all non standard columns
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
        logger.error(f"The following columns could not be retrieved from the annotation file : {columns_not_found}")
        sys.exit(1)

    return columns_indices, columns_names


def check_columns(input_columns, custom_columns):
    """
    Check whether some of the columns wanted by the user are already in the input file

    Args:
        input_columns (str): input data frame columns
        custom_columns (str): annotation columns

    """

    logger = logging.getLogger("logger")

    shared_columns = set(input_columns) & set(custom_columns)

    for shared_column in shared_columns:
        logger.warning(
            f"{shared_column} column exists already in variant file. The one coming from the custom annotation file will be suffixed with _custom"
        )

    return shared_columns


@click.command()
@click.argument("variants", type=click.Path(exists=True))
@click.argument("annotation", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
@click.option(
    "-c", "--columns", default="", type=str, help="Columns to use in the annotation file. Should be comma separated"
)
@click.option(
    "-C", "--columns-file", default="", type=str, help="File listing columns to extract from the annotation file"
)
def cli(variants, annotation, output, columns, columns_file):
    """
    Annotate a TSV variants file (DNM or rates) with information coming from a
    TAB-delimited genome position file indexed with tabix
    (https://www.htslib.org/doc/tabix.html)

    Args:
        variants (str): Variants file (e.g. rates, DNM) starting with the following columns [gene_id,chrom,pos,ref,alt]
        annotation (str): annotation file
        output (str): variants file annotated
        columns (str): columns to use in the annotation file (comma separated)
        columns_file (str): File listing which annotations to extract from dbNSFP

    """

    init_log()

    global is_chr_prefixed_variants, is_chr_prefixed_annotation
    is_chr_prefixed_variants = is_chr_prefixed(variants)
    is_chr_prefixed_annotation = is_chr_prefixed(annotation)

    # Load variants file
    variants_df = load_variants_file(variants)

    # Load indexed annotation file
    annotation_df = pysam.TabixFile(annotation)

    # Get colums
    columns_indices, columns_names = extract_columns_indices(annotation_df, columns, columns_file)

    # If some columns are already found in the input file, we prefix the new ones with "dbnsfp"
    existing_columns = check_columns(variants_df.columns, columns_names)
    columns_names = [f"{col}_custom" if col in existing_columns else col for col in columns_names]

    # Annotate variants file
    annotated_df = annotate(variants_df, annotation_df, columns_indices, columns_names)

    # Export annotated file
    annotated_df.to_csv(output, sep="\t", index=False, na_rep=".")


if __name__ == "__main__":
    cli()
