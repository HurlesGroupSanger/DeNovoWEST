#!/usr/bin/env python
"""Annotate a variant table from a generic tabix-indexed annotation source.

This command is the most flexible annotation entrypoint in the package. It
retrieves one or more user-selected columns from a tabix-indexed TSV-like file
and merges them back onto a variant table using chromosome, position, reference,
and alternate allele.
"""
import gzip
import logging
import sys
from itertools import count, groupby

import click
import pandas as pd
import pysam

from denovowest.utils.io_helpers import as_range, is_chr_prefixed, read_columns_from_file
from denovowest.utils.log import init_log

# TODO : Handle annotation with matching on gene identifier/symbol


def load_variants_file(variants: str) -> pd.DataFrame:
    """Load a variants file (e.g. DNM/rates file).

    Args:
        variants: Path to the variants file.

    Returns:
        DataFrame with variant records.
    """

    return pd.read_table(variants, dtype={"chrom": str, "pos": int, "ref": str, "alt": str})


def annotate(
    variants_df: pd.DataFrame,
    annotation_df: pysam.TabixFile,
    columns_indices: list[int],
    columns_names: list[str],
    chr_prefixed: bool,
) -> pd.DataFrame:
    """Annotate a variant table with columns extracted from a custom resource.

    The annotation source is expected to be tabix-indexed and to use the first
    four columns as ``chrom``, ``pos``, ``ref``, and ``alt``. Matching is done on
    those four fields after chromosome-prefix harmonisation.

    Args:
        variants_df: Input variants to annotate.
        annotation_df: Indexed annotation resource.
        columns_indices: Zero-based column indices to extract from the resource.
        columns_names: Output names for the extracted columns.
        chr_prefixed: Whether chromosome identifiers in the variants file use the "chr" prefix.

    Returns:
        Annotated variant DataFrame.
    """

    list_annotated_df = []
    for gene_id, gene_df in variants_df.groupby("gene_id"):
        annotated_gene_df = annotate_gene(gene_id, gene_df, annotation_df, columns_indices, chr_prefixed)
        list_annotated_df.append(annotated_gene_df)

    # Combine results for all genes and rename columns to their final names
    annotated_df = pd.concat(list_annotated_df)
    annotated_df.rename(dict(zip([str(x) for x in columns_indices], columns_names)), axis=1, inplace=True)

    return annotated_df


def annotate_gene(
    gene_id: str,
    gene_df: pd.DataFrame,
    annotation_df: pysam.TabixFile,
    columns: list[int],
    chr_prefixed: bool,
) -> pd.DataFrame:
    """Annotate variants for a single gene.

    Args:
        gene_id: Gene identifier.
        gene_df: Variants for this gene.
        annotation_df: Tabix-indexed annotation file.
        columns: Indices of annotation columns to retrieve.
        chr_prefixed: Whether chromosome identifiers use the "chr" prefix.

    Returns:
        Annotated variant DataFrame for this gene.
    """

    logger = logging.getLogger("logger")

    gene_annotation_df = retrieve_annotation(gene_df, annotation_df, columns, gene_id, chr_prefixed)

    if not gene_annotation_df.empty:
        gene_annotated_df = gene_df.merge(gene_annotation_df, how="left", on=["chrom", "pos", "ref", "alt"])
    else:
        logger.warning(f"No annotations for gene {gene_id}")
        gene_annotated_df = gene_df

    if gene_df.shape[0] != gene_annotated_df.shape[0]:
        raise ValueError(
            f"Row count changed after merge for {gene_id}: "
            f"{gene_df.shape[0]} -> {gene_annotated_df.shape[0]}"
        )

    return gene_annotated_df


def retrieve_annotation(
    gene_df: pd.DataFrame,
    annotation_df: pysam.TabixFile,
    columns: list[int],
    gene_id: str,
    chr_prefixed: bool,
) -> pd.DataFrame:
    """Retrieve annotations overlapping the genomic span of one gene.

    The function queries the tabix-indexed annotation in contiguous genomic
    blocks built from the observed variant positions, concatenates the matching
    records, normalises chromosome prefixes when needed, and removes duplicate
    variant records before merging back onto the input table.

    Args:
        gene_df: Input variants for one gene.
        annotation_df: Indexed annotation resource.
        columns: Column indices to extract from the annotation resource.
        gene_id: Gene identifier used for logging only.
        chr_prefixed: Whether chromosome identifiers in the variants file use the "chr" prefix.

    Returns:
        DataFrame of annotation records, or an empty DataFrame if none found.
    """

    logger = logging.getLogger("logger")

    gene_chrom = str(gene_df.chrom.values[0]).replace("chr", "")

    # Split each gene in contiguous blocks (i.e. exons) and retrieve annotations
    list_block_df = []
    for _, block in groupby(sorted(set(gene_df["pos"])), key=lambda n, c=count(): n - next(c)):

        # Get the coordinates of the current block
        start, end = as_range(block)

        # Retrieve annotation for the current block
        try:
            block_annotation_df = load_annotation(annotation_df, gene_chrom, start - 1, end, columns)
            list_block_df.append(block_annotation_df)
        except ValueError:
            continue

    if not list_block_df:
        return pd.DataFrame()

    # Concatenate into gene-wide annotation
    gene_annotations_df = pd.concat(list_block_df)

    # Restore the "chr" prefix if the variants file uses it
    if not gene_annotations_df.empty and chr_prefixed:
        gene_annotations_df["chrom"] = "chr" + gene_annotations_df["chrom"]

    # Remove duplicated rows if any
    if not gene_annotations_df.empty:
        n_before = gene_annotations_df.shape[0]
        gene_annotations_df = gene_annotations_df.loc[
            ~gene_annotations_df[["chrom", "pos", "ref", "alt"]].duplicated()
        ]
        if n_before != gene_annotations_df.shape[0]:
            logger.warning(f"Multiple annotations found for {gene_id}; duplicates removed")

    return gene_annotations_df


def load_annotation(
    annotation_file: pysam.TabixFile,
    chrom: str,
    start: int,
    end: int,
    columns: list[int],
) -> pd.DataFrame:
    """Fetch records from a tabix-indexed annotation file for a genomic region.

    Args:
        annotation_file: Tabix-indexed annotation file.
        chrom: Chromosome (without "chr" prefix).
        start: 0-based start position.
        end: End position.
        columns: Column indices to extract from each record.

    Returns:
        DataFrame of annotation records.
    """

    list_records = [_parse(record, columns) for record in annotation_file.fetch(chrom, start, end)]
    records_df = pd.DataFrame(list_records)

    # TODO : handle annotation duplication (especially in the context of gene/transcript id)
    records_df.drop_duplicates(inplace=True)

    return records_df


def _parse(line: str, columns: list[int]) -> dict:
    record = line.split("\t")
    record_dict = {
        "chrom": record[0],
        "pos": int(record[1]),
        "ref": record[2],
        "alt": record[3],
    }
    for column in columns:
        record_dict[str(column)] = record[column]
    return record_dict


def extract_columns_indices(
    annotation_df: pysam.TabixFile,
    columns: str,
    columns_file: str,
) -> tuple[list[int], list[str]]:
    """Extract the indices and names of columns to use from a tabix-indexed annotation file.

    Args:
        annotation_df: Tabix-indexed annotation file.
        columns: Comma-separated column names to extract.
        columns_file: Path to a file listing column names to extract (one per line).

    Returns:
        Tuple of (column_indices, column_names).
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
        columns_to_extract = read_columns_from_file(columns_file)
    # If provided as a string we split it
    elif columns:
        columns_to_extract = columns.split(",")
    # Default: all columns beyond the standard [chrom, pos, ref, alt]
    else:
        columns_to_extract = all_columns[4:]

    # Map names to their positional indices
    columns_indices = []
    columns_names = []
    for idx, name in enumerate(all_columns):
        if name in columns_to_extract:
            columns_indices.append(idx)
            columns_names.append(name)

    columns_not_found = set(columns_to_extract) - set(columns_names)
    if columns_not_found:
        logger.error(f"Columns not found in annotation file: {columns_not_found}")
        sys.exit(1)

    return columns_indices, columns_names


def check_columns(input_columns: list[str], custom_columns: list[str]) -> set[str]:
    """Identify annotation columns that already exist in the variants file.

    Args:
        input_columns: Columns already present in the input variants DataFrame.
        custom_columns: Columns from the annotation file.

    Returns:
        Set of column names present in both.
    """

    logger = logging.getLogger("logger")

    shared_columns = set(input_columns) & set(custom_columns)
    for col in shared_columns:
        logger.warning(
            f"{col} already exists in the variant file; "
            "the annotation column will be suffixed with _custom"
        )
    return shared_columns


@click.command()
@click.argument("variants", type=click.Path(exists=True))
@click.argument("annotation", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
@click.option(
    "-c", "--columns", default="", type=str, help="Columns to extract from the annotation file (comma-separated)"
)
@click.option(
    "-C", "--columns-file", default="", type=str, help="File listing columns to extract from the annotation file"
)
def cli(variants: str, annotation: str, output: str, columns: str, columns_file: str) -> None:
    """Annotate a TSV variants file with columns from a tabix-indexed annotation file.

    The annotation file must be tab-delimited with columns
    [chrom, pos, ref, alt, ...] and indexed with tabix
    (https://www.htslib.org/doc/tabix.html).

    Args:
        variants: Variants file (e.g. rates, DNM) with columns [gene_id, chrom, pos, ref, alt, ...].
        annotation: Tabix-indexed annotation file.
        output: Output annotated variants file.
        columns: Annotation columns to extract (comma-separated).
        columns_file: File listing which annotation columns to extract.
    """

    init_log()

    chr_prefixed_variants = is_chr_prefixed(variants)

    # Load variants file
    variants_df = load_variants_file(variants)

    # Load indexed annotation file
    annotation_df = pysam.TabixFile(annotation)

    # Resolve columns to extract
    columns_indices, columns_names = extract_columns_indices(annotation_df, columns, columns_file)

    # Suffix any column names that clash with the input file
    existing_columns = check_columns(list(variants_df.columns), columns_names)
    columns_names = [f"{col}_custom" if col in existing_columns else col for col in columns_names]

    # Annotate and export
    annotated_df = annotate(variants_df, annotation_df, columns_indices, columns_names, chr_prefixed_variants)
    annotated_df.to_csv(output, sep="\t", index=False, na_rep=".")


if __name__ == "__main__":
    cli()
