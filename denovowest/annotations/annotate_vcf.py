#!/usr/bin/env python
"""Annotate a variant table from a VCF resource.

The VCF annotation command is designed for resources whose annotations are stored
in INFO fields rather than in a TSV-style lookup table. Records are matched by
coordinate and allele, with optional filtering by gene identifier or symbol when
that information is available both in the VCF and in the input data.
"""
import logging
from itertools import count, groupby

import click
import pandas as pd
import pysam

from denovowest.utils.io_helpers import as_range, is_chr_prefixed, load_gff, read_columns_from_file
from denovowest.utils.log import init_log


def load_variants_file(variants: str) -> pd.DataFrame:
    """Load a variants file.

    Args:
        variants: Path to the variants file.

    Returns:
        DataFrame with variant records.
    """

    return pd.read_table(variants, dtype={"chrom": str, "pos": int, "ref": str, "alt": str})


def annotate(
    variants_df: pd.DataFrame,
    annotation_vcf: pysam.VariantFile,
    columns: list[str],
    match_gene: bool,
    gene_mapping: dict[str, str],
    chr_prefixed: bool,
    gene_field: str,
    gene_content: str,
) -> pd.DataFrame:
    """Annotate a variant table with selected INFO fields from a VCF.

    The input table is processed gene by gene to limit the size of each genomic
    query. After retrieval, only the requested INFO fields are appended to the
    existing variant columns.

    Args:
        variants_df: Input variant table.
        annotation_vcf: Indexed VCF resource.
        columns: INFO field names to copy from the VCF.
        match_gene: Whether to filter coordinate matches by gene identity.
        gene_mapping: Mapping from input gene identifiers to symbols when needed.
        chr_prefixed: Whether chromosome identifiers use the "chr" prefix.
        gene_field: INFO field storing the gene identifier in the VCF.
        gene_content: Whether the gene field contains "ensembl_id" or "symbol".

    Returns:
        Annotated variant DataFrame.
    """

    logger = logging.getLogger("logger")

    # We keep track of the columns already in the variant file
    existing_columns = list(variants_df.columns)

    # We loop through each gene in the variant file
    list_annotated_df = []
    for gene_id, gene_df in variants_df.groupby("gene_id"):

        logger.info(f"Annotating {gene_id}")

        # Annotate it with the information from the VCF
        annotated_gene_df = annotate_gene(
            gene_id, gene_df, annotation_vcf, match_gene, gene_mapping, columns,
            chr_prefixed, gene_field, gene_content,
        )

        # Add only the information wanted by the user
        list_annotated_df.append(annotated_gene_df[existing_columns + columns])

    # Combine results for all genes
    return pd.concat(list_annotated_df)


def annotate_gene(
    gene_id: str,
    gene_df: pd.DataFrame,
    annotation_vcf: pysam.VariantFile,
    match_gene: bool,
    gene_mapping: dict[str, str],
    columns: list[str],
    chr_prefixed: bool,
    gene_field: str,
    gene_content: str,
) -> pd.DataFrame:
    """Annotate the variants from a given gene.

    Args:
        gene_id: Gene identifier.
        gene_df: Variant file subsetted to a given gene.
        annotation_vcf: VCF file containing the annotations.
        match_gene: Whether to filter annotations by gene identifier.
        gene_mapping: Maps ENSEMBL gene ID to gene symbol.
        columns: INFO fields to retrieve.
        chr_prefixed: Whether chromosome identifiers use the "chr" prefix.
        gene_field: INFO field storing the gene identifier in the VCF.
        gene_content: Whether the gene field contains "ensembl_id" or "symbol".

    Returns:
        Annotated variant DataFrame for this gene.
    """

    logger = logging.getLogger("logger")

    # We retrieve the annotation in the VCF for this gene
    gene_annotation_df = retrieve_annotation(
        gene_id, gene_df, annotation_vcf, match_gene, gene_mapping,
        chr_prefixed, gene_field, gene_content,
    )

    # We merge it with the variant file
    if not gene_annotation_df.empty:
        gene_annotation_df = gene_annotation_df[["chrom", "pos", "ref", "alt"] + columns]
        gene_annotated_df = gene_df.merge(gene_annotation_df, how="left", on=["chrom", "pos", "ref", "alt"])
        assert gene_df.shape[0] == gene_annotated_df.shape[0]
    else:
        logger.warning(f"No annotations for gene {gene_id}")
        gene_annotated_df = gene_df.copy()
        for column in columns:
            gene_annotated_df[column] = "."

    return gene_annotated_df


def retrieve_annotation(
    gene_id: str,
    gene_df: pd.DataFrame,
    annotation_vcf: pysam.VariantFile,
    match_gene: bool,
    gene_mapping: dict[str, str],
    chr_prefixed: bool,
    gene_field: str,
    gene_content: str,
) -> pd.DataFrame:
    """Retrieve all VCF records overlapping the current gene window.

    The variant positions in ``gene_df`` are collapsed into contiguous genomic
    blocks to reduce the number of random VCF queries. Retrieved records are then
    normalised and filtered by ``format_gene_annotation`` before being merged back
    onto the input variants.

    Args:
        gene_id: Gene identifier used when gene-aware filtering is enabled.
        gene_df: Input variants for the current gene.
        annotation_vcf: Indexed VCF resource.
        match_gene: Whether to filter coordinate matches by gene identity.
        gene_mapping: Mapping from input gene identifiers to symbols when needed.
        chr_prefixed: Whether chromosome identifiers use the "chr" prefix.
        gene_field: INFO field storing the gene identifier in the VCF.
        gene_content: Whether the gene field contains "ensembl_id" or "symbol".

    Returns:
        DataFrame of annotation records, or an empty DataFrame if none found.
    """

    gene_chrom = str(gene_df.chrom.values[0]).replace("chr", "")

    # Split each gene in contiguous blocks (i.e. exons) and retrieve annotations
    list_block_df = []
    for _, block in groupby(sorted(set(gene_df["pos"])), key=lambda n, c=count(): n - next(c)):

        # Get the coordinates of the current block
        start, end = as_range(block)
        # Retrieve annotation for the current block
        try:
            block_df = load_annotation(annotation_vcf, gene_chrom, start - 1, end)
            list_block_df.append(block_df)
        except ValueError:
            continue

    gene_annotations_df = format_gene_annotation(
        gene_id, list_block_df, match_gene, gene_mapping, chr_prefixed, gene_field, gene_content,
    )

    return gene_annotations_df


def format_gene_annotation(
    gene_id: str,
    list_block_df: list[pd.DataFrame],
    match_gene: bool,
    gene_mapping: dict[str, str],
    chr_prefixed: bool,
    gene_field: str,
    gene_content: str,
) -> pd.DataFrame:
    """Format the gene annotation data frame.

    Concatenates per-block annotation DataFrames, applies chromosome prefix
    harmonisation, optionally filters records to those matching the gene
    identifier, and deduplicates by (pos, ref, alt).

    Args:
        gene_id: Gene identifier from the variants file.
        list_block_df: Per-exon-block annotation DataFrames.
        match_gene: Whether to filter annotations by gene identifier.
        gene_mapping: Maps ENSEMBL gene ID to gene symbol.
        chr_prefixed: Whether chromosome identifiers use the "chr" prefix.
        gene_field: INFO field storing the gene identifier in the VCF.
        gene_content: Whether the gene field contains "ensembl_id" or "symbol".

    Returns:
        Formatted annotation DataFrame, or an empty DataFrame if none found.
    """

    if not list_block_df:
        return pd.DataFrame()

    # Concatenate into gene-wide annotation
    gene_annotations_df = pd.concat(list_block_df)

    if gene_annotations_df.empty:
        return pd.DataFrame()

    # Add the chr prefix if found in the variants file
    if chr_prefixed:
        gene_annotations_df["chrom"] = "chr" + gene_annotations_df["chrom"]

    # Keep only annotations that match the gene identifier
    # TODO : there are a few caveats. It will work only if the variants file contains an ENSEMBL identifier
    # as gene_id and if the VCF contains either ensembl id or gene symbol
    if match_gene:
        # Strip the version suffix from gene_id before matching (handles both "ENSG1" and "ENSG1.7")
        bare_gene_id = gene_id.split(".")[0]
        if gene_content == "ensembl_id":
            gene_annotations_df = gene_annotations_df.loc[gene_annotations_df[gene_field] == bare_gene_id]
        else:
            gene_symbol = gene_mapping[bare_gene_id]
            gene_annotations_df = gene_annotations_df.loc[gene_annotations_df[gene_field] == gene_symbol]

    # TODO In some cases we have one record per transcript, for now we just pick the first record but ideally we could match on transcript
    gene_annotations_df = gene_annotations_df.loc[~gene_annotations_df[["pos", "ref", "alt"]].duplicated()]

    return gene_annotations_df


def load_annotation(annotation_vcf: pysam.VariantFile, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Fetch records from a VCF annotation file for a genomic region.

    Args:
        annotation_vcf: VCF file containing variant annotations.
        chrom: Chromosome identifier (without "chr" prefix).
        start: 0-based start position.
        end: End position.

    Returns:
        DataFrame of annotation records.
    """

    list_records = []
    for record in annotation_vcf.fetch(chrom, start, end):
        list_records.append(_parse(record))

    records_df = pd.DataFrame(list_records)

    # TODO : handle annotation duplication (especially in the context of gene/transcript id)
    # records_df.drop_duplicates(inplace=True)

    return records_df


def _parse(record) -> dict:
    """Turn a VCF record into a dictionary.

    Args:
        record: pysam VCF record.

    Returns:
        Dictionary with variant coordinates and INFO fields.
    """

    record_dict = {"chrom": record.chrom, "pos": record.pos, "ref": record.ref, "alt": record.alts[0]}
    record_dict = record_dict | dict(record.info)
    return record_dict


def detect_gene_field_vcf(annotation: pysam.VariantFile) -> tuple[str, str]:
    """Detect which INFO field in the VCF stores the gene identifier or symbol.

    Args:
        annotation: VCF annotation file.

    Returns:
        Tuple of (gene_field, gene_content) where gene_content is "ensembl_id"
        or "symbol", or ("", "") if no usable gene field is found.
    """

    logger = logging.getLogger("logger")

    record = next(iter(annotation))
    record_info = dict(record.info)

    gene_fields = list(set(["gene", "GENE"]) & set(record_info.keys()))

    if not gene_fields:
        logger.warning("Could not find a gene field in the VCF INFO fields. Matching variants on coordinates only")
        return "", ""

    if len(gene_fields) > 1:
        logger.warning("Several fields can contain the gene information. Matching variants on coordinates only")
        return "", ""

    gene_field = gene_fields[0]
    gene_content = "ensembl_id" if record_info[gene_field].startswith("ENSG") else "symbol"

    return gene_field, gene_content


def build_gene_mapping_from_gff(gff_db) -> dict[str, str]:
    """Build mapping between ENSEMBL gene identifier and gene symbol.

    #TODO : handle cases with multiple name/ids

    Args:
        gff_db: GFF database.

    Returns:
        Dictionary mapping ENSEMBL gene ID to gene symbol.
    """

    gene_mapping = {}
    for gene in gff_db.all_features(featuretype="gene"):

        gene_id = gene.attributes["gene_id"][0].split(".")[0]

        if "gene_name" in dict(gene.attributes).keys():
            gene_name = gene.attributes["gene_name"][0]
        elif "Name" in dict(gene.attributes).keys():
            gene_name = gene.attributes["Name"][0]
        else:
            gene_name = ""

        gene_mapping[gene_id] = gene_name

    return gene_mapping


def check_columns(
    variants_df: pd.DataFrame,
    annotation_vcf: pysam.VariantFile,
    columns: str,
    columns_file: str,
) -> list[str]:
    """Determine which INFO fields to retrieve from the VCF annotation file.

    Args:
        variants_df: Variants file.
        annotation_vcf: VCF file containing the annotations.
        columns: Comma-separated INFO field names to retrieve.
        columns_file: File listing INFO field names to retrieve (one per line).

    Returns:
        List of INFO field names to retrieve.
    """

    logger = logging.getLogger("logger")

    # If the user provided the columns to extract in a separate file we read it
    if columns_file:
        annotation_columns = read_columns_from_file(columns_file)
    # If provided as a string we split it
    elif columns:
        annotation_columns = columns.split(",")
    # Default: all INFO fields from the first record
    else:
        record = next(iter(annotation_vcf))
        annotation_columns = list(dict(record.info).keys())

    shared_columns = set(variants_df.columns) & set(annotation_columns)
    for col in shared_columns:
        logger.warning(f"{col} already exists in the variant file")

    return annotation_columns


def build_gene_mapping_from_column(variants_df: pd.DataFrame) -> dict[str, str]:
    """Build a dictionary that maps ENSEMBL gene identifiers to gene symbols.

    Args:
        variants_df: Variants file containing gene_id and symbol columns.

    Returns:
        Dictionary mapping ENSEMBL gene ID to gene symbol.
    """

    df = variants_df[["gene_id", "symbol"]].drop_duplicates(subset=["gene_id"])
    return dict(zip(df.gene_id, df.symbol))


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
@click.option("--match-gene", is_flag=True, help="Match the annotation based on gene identifier or symbol")
@click.option(
    "--gff",
    type=str,
    help="GFF file used to build the rates file, and that contains the ENSEMBL identifier / gene symbol matching",
)
def cli(variants: str, annotation: str, output: str, columns: str, columns_file: str, match_gene: bool, gff: str) -> None:
    """Annotate a variants file (e.g. rates file) with informations stored in a VCF.

    CAVEAT: For now the script works only with variant files using ENSEMBL gene identifiers,
    and with VCFs such as the popeve and Cosmic ones. There is no guarantee that it will work
    with other data.

    Args:
        variants: Variant file.
        annotation: Annotation to add to the variant file in VCF format.
        output: Annotated variant file.
        columns: Columns to retrieve from the VCF (comma-separated).
        columns_file: File listing which annotations to extract from the VCF.
        match_gene: Whether to match annotation using coordinates only or gene id.
        gff: GFF file used to build the rates file; contains matching between gene
            identifier and gene symbol that can be used if the annotation contains
            only gene symbol.
    """

    init_log()

    # Handle differences in chromosome representations
    chr_prefixed_variants = is_chr_prefixed(variants)

    # Load variants file
    variants_df = load_variants_file(variants)

    # Load indexed annotation file
    annotation_vcf = pysam.VariantFile(annotation)

    # See if we can retrieve gene information in the VCF
    gene_field, gene_content = detect_gene_field_vcf(annotation_vcf)
    if not gene_field:
        match_gene = False

    # Some annotations rely on gene symbol rather than identifiers; build a mapping dictionary
    # either from the GFF file that was used to generate the rates file
    if gff:
        gff_db = load_gff(gff)
        gene_mapping = build_gene_mapping_from_gff(gff_db)
    # or from a symbol column if found in the file
    elif "symbol" in variants_df.columns:
        gene_mapping = build_gene_mapping_from_column(variants_df)
    else:
        gene_mapping = {}

    # Annotate variants file
    annotation_columns = check_columns(variants_df, annotation_vcf, columns, columns_file)
    annotated_df = annotate(
        variants_df, annotation_vcf, annotation_columns, match_gene, gene_mapping,
        chr_prefixed_variants, gene_field, gene_content,
    )

    # Export annotated file
    annotated_df.to_csv(output, sep="\t", index=False, na_rep=".")


if __name__ == "__main__":
    cli()
