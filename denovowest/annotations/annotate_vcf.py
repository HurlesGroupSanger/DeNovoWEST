#!/usr/bin/env python
import click
import pysam
import pandas as pd
from itertools import groupby, count
import gffutils
import os
import logging
import utils


def is_chr_prefixed(df_path):
    """
    Check whether the chromosomes are prefixed with "chr" in the input tables

    Args:
        df_path (str): path to a variant or annotation file

    """

    df = pd.read_csv(df_path, sep="\t", comment="#", nrows=10000)

    # We first check for a chrom column
    for key in ["chrom", "chr", "#chr", "#chrom"]:
        if key in df.columns:
            is_chr_prefixed = str(df[key].iloc[0]).startswith("chr")
            break

    # If there is no chrom column we try with the first one
    # TODO : not optimal
    if "is_chr_prefixed" not in locals():
        is_chr_prefixed = str(df.iloc[0, 0]).startswith("chr")

    return is_chr_prefixed


def load_variants_file(variants):
    """
    Load variants file

    Args:
        variants (str): variants file
    """

    variants_df = pd.read_table(variants, dtype={"chrom": str, "pos": int, "ref": str, "alt": str})

    return variants_df


def annotate(variants_df, annotation_vcf, columns, match_gene, gene_mapping):
    """
    Annotate the variants file with the informations from VCF

    Args:
        variants_df (pd.DataFrame): variants file
        annotation_vcf (VariantFile): VCF file containing the annotations
        columns (list): info fields to retrieve from the VCF

    """

    logger = logging.getLogger("logger")

    # We keep track of the columns already in the variant file
    existing_columns = list(variants_df.columns)

    # We loop through each gene in the variant file
    list_annotated_df = list()
    for gene_id, gene_df in variants_df.groupby("gene_id"):

        logger.info(f"Annotating {gene_id}")

        # Annotate it with the informations from the VCF
        annotated_gene_df = annotate_gene(gene_id, gene_df, annotation_vcf, match_gene, gene_mapping, columns)

        # Add only the informations wanted by the user
        list_annotated_df.append(annotated_gene_df[existing_columns + columns])

    # Combine results for all genes
    annotated_df = pd.concat(list_annotated_df)

    return annotated_df


def annotate_gene(gene_id, gene_df, annotation_vcf, match_gene, gene_mapping, columns):
    """
    Annotate the variants from a given gene

    Args:
        gene_id (str): gene identifier
        gene_df (pd.DataFrame): variant file subsetted to a given gene
        annotation_vcf (VariantFile): VCF file containing the annotations

    """

    logger = logging.getLogger("logger")

    # We retrieve the annotation in the VCF for this gene
    gene_annotation_df = retrieve_annotation(gene_id, gene_df, annotation_vcf, match_gene, gene_mapping)

    # We merge it with the variant file
    if not gene_annotation_df.empty:
        gene_annotated_df = gene_df.merge(gene_annotation_df, how="left", on=["chrom", "pos", "ref", "alt"])
        assert gene_df.shape[0] == gene_annotated_df.shape[0]
    else:
        logger.warning(f"No annotations for gene {gene_id}")

        gene_annotated_df = gene_df
        for column in columns:
            gene_annotated_df[column] = "."

    return gene_annotated_df


def retrieve_annotation(gene_id, gene_df, annotation_vcf, match_gene, gene_mapping):
    """
    Retrieve annotation from the VCF file

    Args:
        gene_df (pd.DataFrame): variant file subsetted to a given gene
        annotation_vcf (VariantFile): VCF file containing the annotations

    """

    gene_chrom = str(gene_df.chrom.values[0]).replace("chr", "")

    # Split each gene in contiguous block (i.e. exons) and retrieve annotations
    list_block_df = list()
    for _, block in groupby(sorted(set(gene_df["pos"])), key=lambda n, c=count(): n - next(c)):

        # Get the coordinates of the current block
        start, end = as_range(block)
        # Retrieve annotation for the current block
        try:
            block_df = load_annotation(annotation_vcf, gene_chrom, start - 1, end)
            list_block_df.append(block_df)
        except ValueError as e:
            continue

    gene_annotations_df = format_gene_annotation(gene_id, list_block_df, match_gene, gene_mapping)

    return gene_annotations_df


def format_gene_annotation(gene_id, list_block_df, match_gene, gene_mapping):
    """
    Format the gene annotation data frame

    Args:
        gene_id (_type_): _description_
        list_block_df (_type_): _description_
        match_gene (_type_): _description_
        gff (_type_): _description_
    """

    if len(list_block_df) == 0:
        return pd.DataFrame()

    # Concatenate into gene wide annotation
    gene_annotations_df = pd.concat(list_block_df)

    if gene_annotations_df.empty:
        return pd.DataFrame()

    # Add the chr prefix if found in the variants file
    if is_chr_prefixed_variants:
        gene_annotations_df["chrom"] = ["chr" + x for x in gene_annotations_df["chrom"]]

    # Keep only annotations that match the gene identifier
    # TODO : there are a few caveats. It will work only if the variants file contains an ENSEMBL identifier
    # as gene_id and if the VCF contains either ensembl id or gene symbol
    if match_gene:

        if gene_content == "ensembl_id":
            # TODO : handle ensembl gene version
            gene_annotations_df = gene_annotations_df.loc[gene_annotations_df[gene_field] == gene_id]
        else:
            gene_symbol = gene_mapping[gene_id]
            gene_annotations_df = gene_annotations_df.loc[gene_annotations_df[gene_field] == gene_symbol]

    # TODO In some cases we have one record per transcript, for now we just pick the first record but ideally we could match on transcript
    gene_annotations_df = gene_annotations_df.loc[~gene_annotations_df[["pos", "ref", "alt"]].duplicated()]

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


def load_annotation(annotation_vcf, chrom, start, end):
    """
    Access to specific region in annotation file
    Args:
        annotation_vcf (pysam.VariantFile): VCF file containing variant annotations
        chrom (str): chromosome identifier
        start (int): beginning position of the region
        end (int): terminating position of the region
    """

    list_records = list()
    for record in annotation_vcf.fetch(chrom, start, end):
        list_records.append(parse(record))

    records_df = pd.DataFrame(list_records)

    # TODO : handle annotation duplication (especially in the context of gene/transcript id)
    # records_df.drop_duplicates(inplace=True)

    return records_df


def parse(record):
    """
    Turns the VCF record into a dictionnary

    Args:
        record (VCFRecord): VCF record

    """

    record_dict = {"chrom": record.chrom, "pos": record.pos, "ref": record.ref, "alt": record.alts[0]}
    record_dict = record_dict | dict(record.info)
    return record_dict


def detect_gene_field_vcf(annotation):
    """
    Because we may want to match variants on gene identifiers, we need to retrieve
    which field in the VCF stores the gene identifier/symbol.

    Args:
        annotation (VariantFile): variant file
    """

    logger = logging.getLogger("logger")

    for i, record in enumerate(annotation):
        break

    record_info = dict(record.info)

    gene_fields = list(set(["gene", "GENE"]) & set(record_info.keys()))

    if gene_fields:

        if len(gene_fields) > 1:
            logger.warning("Several fields can contain the gene information. Matching variants on coordinates only")
        else:
            gene_field = gene_fields[0]
            if record_info[gene_field].startswith("ENSG"):
                gene_content = "ensembl_id"
            else:
                gene_content = "symbol"
    else:
        gene_field = ""
        gene_content = ""
        logger.warning("Could not find a gene field in the VCF INFO fields. Matching variants on coordinates only")

    return gene_field, gene_content


def load_gff(gff_file, gff_db_path):
    """Create the gff database used by gffutils.

    Args:
        gff_file (str): Path to a GFF or a gffutils database file.
        gff_db_path (str): Path to the newly created gff database.

    Returns:
        gffutils.db: GFF database.

    """

    logger = logging.getLogger("logger")

    # gffutils db input
    if gff_file.endswith(".db"):
        logger.info(f"Loading gffutils database : {gff_file}")
        gff_db = gffutils.FeatureDB(gff_file)
    # GFF input
    else:
        try:
            os.remove(gff_db_path)
            logger.info(f"Removed old gffutil database : {gff_db_path}")
        except OSError:
            pass

        logger.info(f"Creating gffutil database : {gff_db_path}")
        gff_db = gffutils.create_db(gff_file, gff_db_path, merge_strategy="create_unique")

    return gff_db


def build_gene_mapping_from_gff(gff_db):
    """
    Build mapping between gene identifier (ENSEMBL) and symbol.
    #TODO : handle cases with multiple name/ids

    Args:
        gff_db (gffutils.db): GFF database.

    """

    gene_mapping = dict()
    for gene in gff_db.all_features(featuretype="gene"):

        gene_id = gene.attributes["gene_id"][0]
        gene_name = gene.attributes["gene_name"][0]

        gene_mapping[gene_id] = gene_name

    return gene_mapping


def check_columns(variants_df, annotation_vcf, columns, columns_file):
    """
    Check which column to retrieve from the VCF annotation file

    Args:
        variants_df (pd.DataFrame): variants file
        annotation_vcf (VariantFile): VCF file containing the annotations
        columns (list): list of annotation to retrieve

    """

    logger = logging.getLogger("logger")

    # If the user provided the columns to extract in a separate file we read it
    if columns_file:
        annotation_columns = utils.read_columns_from_file(columns_file)
    # If he provided them as a string we split it
    elif columns:
        annotation_columns = columns.split(",")
    # If he did not provide any column we retrieve all the fields in INFO
    else:
        for i, record in enumerate(annotation_vcf):
            break
        annotation_columns = list(dict(record.info).keys())

    variants_columns = variants_df.columns

    shared_columns = set(variants_columns) & set(annotation_columns)

    for shared_column in shared_columns:
        logger.warning(f"{shared_column} column exists already in variant file")

    return annotation_columns


def build_gene_mapping_from_column(variants_df):
    """
    Build a dictionnary that maps gene ENSEMBL identifiers to gene symbol

    Args:
        variants_df (pd.DataFrame): variants file
    """

    df = variants_df[["gene_id", "symbol"]]
    df = df.loc[~df[["gene_id"]].duplicated()]

    gene_mapping = dict(zip(df.gene_id, df.symbol))

    return gene_mapping


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
def cli(variants, annotation, output, columns, columns_file, match_gene, gff):
    """
    Annotate a variants file (e.g. rates file) with informations stored in a VCF.
    CAVEAT : For now the script works only with variant file using ENSEMBL gene identifiers, and with VCF
    such as the popeve and Cosmic one. There is no guarantee that it will work with other data.

    Args:
        variants (str): variant file
        annotation (str): annotation to add to the variant file in VCF format
        output (str): annotated variant file
        columns (str): columns to retrieve from the VCF. Should be comma separated
        columns_file (str): File listing which annotations to extract from dbNSFP
        match_gene (boolean): whether to match annotation using coordinates only or gene id
        gff (str): gff file used to build the rates file, contains matching between gene identifier and gene symbol that
        can be used if the annotation contains only gene symbol
    """

    utils.init_log()

    # Handle differences in chromosome representations
    global is_chr_prefixed_variants, is_chr_prefixed_annotation
    is_chr_prefixed_variants = is_chr_prefixed(variants)
    is_chr_prefixed_annotation = is_chr_prefixed(annotation)

    # Load variants file
    variants_df = load_variants_file(variants)

    # Load indexed annotation file
    annotation_vcf = pysam.VariantFile(annotation)

    # See if we can retrieve gene informations in the VCF
    global gene_field, gene_content
    gene_field, gene_content = detect_gene_field_vcf(annotation_vcf)
    if gene_field == "":
        match_gene = False

    # Some annotations rely on gene symbol rather than identifiers, we therefore build a mapping dictionnary
    # either from the gff file that was used to generate the rates file
    if gff:
        gff_db = load_gff(gff, "gff.db")
        gene_mapping = build_gene_mapping_from_gff(gff_db)
    # or from a symbol column if found in the file
    elif "symbol" in variants_df.columns:
        gene_mapping = build_gene_mapping_from_column(variants_df)
    else:
        gene_mapping = dict()

    # Annotate variants file
    annotation_columns = check_columns(variants_df, annotation_vcf, columns, columns_file)
    annotated_df = annotate(variants_df, annotation_vcf, annotation_columns, match_gene, gene_mapping)

    # Export annotated file
    annotated_df.to_csv(output, sep="\t", index=False, na_rep=".")


if __name__ == "__main__":
    cli()
