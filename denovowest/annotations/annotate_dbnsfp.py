#!/usr/bin/env python
import pandas as pd
import click
import pysam
import sys
import logging

from itertools import groupby, count
from denovowest.utils.log import init_log
from denovowest.utils.io_helpers import read_columns_from_file, load_gff


# TODO : Better handle the ENSG versions
# TODO : See if we could improve the selection when multiple records exist for a same variant (@remove_duplicates)


def load_dbnsfp(
    dbnsfp_file, chrom, start, end, columns_indices, gene_id, gff_db=None, ensembl_gene_id_map_version=None
):
    """
    Access to specific region in an dbNSFP indexed file
    Args:
        annotation_file (pysam.TabixFile): Annotation file containing CEP scores
        chrom (str): chromosome identifier
        start (int): beginning position of the region
        end (int): terminating position of the region
        columns_indices(list) : indices of the dbNSFP columns to retrieve
        gene_id(str) : gene identifier used in DNM and rates file
        gff_db (gffutils.FeatureDB): gffutils database
        ensembl_gene_id_map_version (dict) : maps ENSG with version to ENSG without version (e.g. {ENSG00000010404 : ENSG00000010404.1})
    """

    # Retrieve the transcripts associated to the genes from the GFF file
    if gff_db:
        transcript_ids = get_transcripts(gene_id, gff_db, ensembl_gene_id_map_version)
    else:
        transcript_ids = list()

    list_annotated_records = list()
    for record in dbnsfp_file.fetch(chrom, start, end):
        record_dict = parse(record, columns_indices, gene_id, transcript_ids)
        if record_dict:
            list_annotated_records.append(record_dict)

    return pd.DataFrame(list_annotated_records)


def get_transcripts(gene_id, gff_db, ensembl_gene_id_map_version):
    """
    Return transcripts identifiers of the transcripts found in the GFF file for a given gene

    Args:
        gene_id (str): gene identifier
        gff_db (gffutils.FeatureDB): gffutils database
        ensembl_gene_id_map_version (dict) : maps ENSG with version to ENSG without version (e.g. {ENSG00000010404 : ENSG00000010404.1})

    Returns:
        list: list of transcripts identifiers
    """

    if "." in gene_id:
        gene = gff_db[gene_id]
    else:
        gene = gff_db[ensembl_gene_id_map_version[gene_id]]

    list_transcript_ids = list()
    for transcript in gff_db.children(gene, featuretype="transcript", order_by="start"):
        list_transcript_ids += [x.split(".")[0] for x in transcript["transcript_id"]]

    return list_transcript_ids


def parse(record, columns_indices, gene_id, transcript_ids_gff=list()):
    """
    Parse a dbNSFP record (line) and extract the data of interest

    Args:
        record (str): line in dbNSFP
        columns_indices (list): indices of columns to retrieve from the dbNSFP records
        gene_id(str) : ENSEMBL gene identifier from the input file
        transcript_ids_gff (list) : list of {gene_id} transcripts found in the GFF file (if provided)

    Returns:
        dict: informations about a mutation (pos, CEP scores)
    """

    record = record.split("\t")

    record_dict = dict()
    record_dict["chrom"] = record[0]
    record_dict["pos"] = int(record[1])
    record_dict["ref"] = record[2]
    record_dict["alt"] = record[3]

    # TODO : should retrieve the index rather than relying ona position that might change with versions of dbNSDP
    ensembl_gene_ids_dbnsfp = record[13].split(";")
    transcript_ids_dbnsfp = record[14].split(";")

    # Check the intersection between the transcripts from the user gff and dbNSFP
    if len(set(transcript_ids_gff) & set(transcript_ids_dbnsfp)) != 0:
        transcript_in_dbnsfp = True
    else:
        transcript_in_dbnsfp = False
    record_dict["transcript_in_dbnsfp"] = transcript_in_dbnsfp

    # If the record is associated to an overlapping gene, skip it
    if gene_id.split(".")[0] not in ensembl_gene_ids_dbnsfp:
        return dict()

    # Retrieve the scores
    for column_index in columns_indices:
        record_dict[str(column_index)] = record[column_index]

    # If a GFF has been provided, retrieve the maximum transcript-based score
    if transcript_ids_gff:

        for column_index in columns_indices:
            list_scores = record[column_index].split(";")

            score = "."

            # If there is only one score, at this point we assume it is a per-gene score
            if len(list_scores) == 1:
                score = list_scores[0]
                record_dict[str(column_index)] = score
                continue

            # We first try to retrieve the transcript based score
            if transcript_in_dbnsfp:
                score = get_transcript_based_score(list_scores, transcript_ids_gff, transcript_ids_dbnsfp)

            # Some MANE transcripts are missing in dbNSFP, therefore if no score associated to a transcript could be retrieved
            # we switch to gene mode
            if score == ".":
                score = get_gene_based_score(list_scores, gene_id.split(".")[0], ensembl_gene_ids_dbnsfp)

            record_dict[str(column_index)] = score

    return record_dict


def get_transcript_based_score(list_scores, transcript_ids_gff, transcript_ids_dbnsfp):
    """
    Return the maximum transcript-based score among transcripts found in the GFF

    Args:
        list_scores (list): list of transcript-based scores for a given metric
        transcript_ids_gff(list) : list of transcript identifiers found in the GFF
        transcript_ids_dbnsfp(list) : list of transcript identifiers in dbNSFP for the given variant
    """

    list_idx = list()
    for idx, transcript_id in enumerate(transcript_ids_dbnsfp):
        if transcript_id in transcript_ids_gff:
            list_idx.append(idx)

    # Not all dbNSFP columns contain scores, and they do not necessarily match the number of values in gene/transcript columns
    try:
        max_score = max_value([list_scores[i] for i in list_idx])
    except IndexError:
        # In that case we just return the whole value
        max_score = ";".join(list_scores)

    return max_score


def get_gene_based_score(list_scores, gene_id, gene_ids_dbnsfp):
    """
    Return the maximum gene-based score among all transcripts associated to this gene

    Args:
        list_scores (list): list of transcript-based scores for a given metric
        transcript_ids_gff(list) : list of transcript identifiers found in the GFF
        transcript_ids_dbnsfp(list) : list of transcript identifiers in dbNSFP for the given variant
    """

    list_idx = list()
    for idx, dbnsfp_gene_id in enumerate(gene_ids_dbnsfp):
        if dbnsfp_gene_id == gene_id:
            list_idx.append(idx)

    # Not all dbNSFP columns contain scores, and they do not necessarily match the number of values in gene/transcript columns
    try:
        max_score = max_value([list_scores[i] for i in list_idx])
    except IndexError:
        # In that case we just return the whole value
        max_score = ";".join(list_scores)

    return max_score


def max_value(list_scores):
    """
    Return the maximum value in a list that can contains "."

    Args:
        list_scores (list): list of transcript-based scores for a given metric

    Returns:
        float: maximum score in the list
    """
    try:
        max_score = max([float(x) for x in list_scores if x != "."])
    except ValueError:

        # Not all dbNSFP columns contain scores, and they do not necessarily match the number of values in gene/transcript columns
        if set(list_scores) == {"."}:
            # If there is no score we return "."
            max_score = "."
        else:
            # If it was not a score column we return the whole value
            max_score = ";".join(list_scores)

    return max_score


def as_range(region):
    """
    Returns genomic region boundaries

    Args:
        region (list): genomic region positions

    Returns:
        tuple: beginning and terminating position of the genomic region
    """

    return region[0], region[-1]


def extract_columns_from_dbNFP(dbnsfp, columns, columns_file):
    """
    Extract the indices of the column to use in dbnsfp

    Args:
        dbnsfp (pysam.TabixFile): dbNSFP indexed file
        columns(str) : Annotations to extract form dbNSFP (comma separated)
        columns_file (str): File listing which annotations to extract from dbNSFP
    """

    logger = logging.getLogger("logger")

    # If the user provided the columns to extract in a separate file we read it
    if columns_file:
        columns_to_extract = read_columns_from_file(columns_file)
    # If he provided them as a string we split it
    elif columns:
        columns_to_extract = columns.split(",")
    # If he did not provide any column we retrieve all columns in dbNSFP
    else:
        columns_to_extract = dbnsfp.header[0].split("\t")[4:]

    # Retrieve dbNSFP columns
    dbnsfp_columns = dbnsfp.header[0].split("\t")

    # Find the index of all columns asked by the user
    list_indices = list()
    found_column = list()
    for idx, name in enumerate(dbnsfp_columns):
        if name in columns_to_extract:
            list_indices.append(idx)
            found_column.append(name)

    # If a column can not be retrieved from dbNSFP we stop here
    columns_not_found = set(columns_to_extract) - set(found_column)
    if columns_not_found:
        logger.error(f"The following columns could not be retrieved from dbNFSP : {columns_not_found}")
        sys.exit(1)

    return list_indices, found_column


def check_columns(input_columns, dbnsfp_columns):
    """
    Check whether some of the columns wanted by the user are already in the input file

    Args:
        input_columns (str): input data frame columns
        dbnsfp_columns (str): annotation columns

    """

    logger = logging.getLogger("logger")

    shared_columns = set(input_columns) & set(dbnsfp_columns)

    for shared_column in shared_columns:
        logger.warning(
            f"{shared_column} column exists already in variant file. The one coming from dbNSFP will be suffixed with _dbnsfp"
        )

    return shared_columns


def load_rates_file(rates, output):
    """
    Load rates file

    Args:
        rates (str): path to a variant file in TSV format
        output (str): output file

    """

    logger = logging.getLogger("logger")

    # Load rates file
    rates_df = pd.read_table(rates, dtype={"chrom": str, "pos": int, "ref": str, "alt": str})

    # Edge case where when splitting the processes we end up with an empty input file
    if rates_df.empty:
        logger.info("Rates file is empty")
        rates_df["raw"] = None
        rates_df["score"] = None
        rates_df.to_csv(output, sep="\t", index=False)
        sys.exit(0)

    # Depending on the gff, chromosome can be defined as "chrN" or just "N"
    if str(rates_df.iloc[0].chrom).startswith("chr"):
        add_chr = True
    else:
        add_chr = False

    return rates_df, add_chr


def remove_duplicates(block_dbnsfp):
    """
    Looking at dbNSFP it looks like we can have two records for a same variant if it leads to
    a different amino acid.
    For now we prioritise the record that matches the transcript in the input GFF (if any) and,
    if several, the record with less missing annotations.


    Args:
        block_dbnsfp (pd.DataFrame): subset of the annotated data frame
    """

    indices_to_remove = list()
    for pos, min_df in block_dbnsfp.loc[block_dbnsfp[["chrom", "pos", "ref", "alt"]].duplicated(keep=False)].groupby(
        ["pos", "alt"]
    ):

        indices = min_df.index
        if min_df["transcript_in_dbnsfp"].sum() != 0:
            min_df = min_df.loc[min_df["transcript_in_dbnsfp"] == True].copy()

        # Keep the record with the most values annotated
        min_df["nb_missing_values"] = (min_df != ".").sum(axis=1)
        keep_idx = min_df.sort_values("nb_missing_values", ascending=False).index[0]

        indices_to_remove = indices_to_remove + list(set(indices) - set([keep_idx]))

    block_dbnsfp = block_dbnsfp.drop(indices_to_remove)

    return block_dbnsfp


def extract_ensembl_gene_id_without_version(gff_db):
    """
    It might happen that the user input file contains ENSEMBL gene ids without version number
    when the GFF file does contain them. In that case we need to map the two together.

    Args:
        gff_db (gffutils.FeatureDB): gffutils database

    Returns:
        dict: maps ENSG with version to ENSG without version (e.g. {ENSG00000010404 : ENSG00000010404.1})
    """

    ensembl_gene_id_map_version = dict()
    for gene in gff_db.features_of_type("gene", order_by="start"):
        ensembl_gene_id_map_version[gene.id.split(".")[0]] = gene.id

    return ensembl_gene_id_map_version


@click.command()
@click.argument("rates_dnm")
@click.argument("dbnsfp")
@click.argument("output")
@click.option("-c", "--columns", default="")
@click.option("-C", "--columns-file", default="")
@click.option("--gff", default="")
def annotate_dbnsfp(rates_dnm, dbnsfp, output, columns, columns_file, gff):
    """
    Annotatate a rates/DNM file with CEP scores from dbNSFP.
    When multiple scores exist for a variant (e.g. multiple transcripts, overlapping genes),
    the maximum score among transcripts found in the GFF are retrieved.
    When no GFF are provided, the maximum score for the gene associated to the variant is retrieved.

    Args:
        rates_dnm (str): Variants file (e.g. rates, DNM) starting with the following columns [gene_id,chrom,pos,ref,alt]
        dbnsfp (str): dbNSFP genome wide file
        output (str): Output file (annotated dataframe)
        columns(str) : Annotations to extract form dbNSFP (comma separated)
        columns_file (str): File listing which annotations to extract from dbNSFP
        gff(str) : Gff file or gffutils database
    """

    init_log()
    logger = logging.getLogger("logger")

    # Load rates/DNM file
    df, add_chr = load_rates_file(rates_dnm, output)
    rates_df_columns = list(df.columns)

    # Load dbnsfp file
    dbnsfp_df = pysam.TabixFile(dbnsfp, encoding="utf-8")

    # Load gff file/DB (to match annotation based on user selected transcripts)
    if gff:
        gff_db = load_gff(gff)
        ensembl_gene_id_map_version = extract_ensembl_gene_id_without_version(gff_db)
    else:
        gff_db = None

    # Retrieve columns to extract from dbnsfp_df
    dbnsfp_columns_indices, dbnsfp_columns_names = extract_columns_from_dbNFP(dbnsfp_df, columns, columns_file)

    # If some columns are already found in the input file, we prefix the new ones with "dbnsfp"
    existing_columns = check_columns(rates_df_columns, dbnsfp_columns_names)
    dbnsfp_columns_names = [f"{col}_dbnsfp" if col in existing_columns else col for col in dbnsfp_columns_names]

    # For each gene
    list_merged_df = list()
    for gene_id, gene_rates_df in df.groupby("gene_id"):

        logger.info(f"Annotating {gene_id}")

        chrom = str(gene_rates_df.chrom.values[0]).replace("chr", "")

        # Split each gene in contiguous block (i.e. exons) and load dbNSFP scores
        list_block_df = list()
        for _, block in groupby(sorted(set(gene_rates_df["pos"])), key=lambda n, c=count(): n - next(c)):
            start, end = as_range(list(block))
            try:
                block_dbnsfp = load_dbnsfp(
                    dbnsfp_df,
                    chrom,
                    start - 1,
                    end,
                    dbnsfp_columns_indices,
                    gene_id,
                    gff_db,
                    ensembl_gene_id_map_version,
                )
                if not block_dbnsfp.empty:
                    block_dbnsfp = remove_duplicates(block_dbnsfp)
                    list_block_df.append(block_dbnsfp)
            except ValueError:
                continue

        # Merge annotations from several blocks together
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
    merged_df.columns = rates_df_columns + ["transcript_in_dbnsfp"] + dbnsfp_columns_names
    merged_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    merged_df = annotate_dbnsfp()
