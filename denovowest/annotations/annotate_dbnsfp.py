#!/usr/bin/env python
import pandas as pd
import click
import pysam
import sys
import gffutils
import os
from itertools import groupby, count


def load_dbnsfp(dbnsfp_file, chrom, start, end, columns_indices, gene_id, gff_db=None):
    """
    Access to specific region in an dbNSFP indexed file
    Args:
        annotation_file (pysam.TabixFile): Annotation file containing CEP scores
        chrom (str): chromosome identifier
        start (int): beginning position of the region
        end (int): terminating position of the region
        columns_indices(list) : indices of the dbNSFP columns to retrieve
        gene_id(str) : gene identifier used in DNM and rates file
        gff_db()
    """

    # Retrieve the transcripts associated to the genes from the GFF file
    if gff_db:
        transcript_ids = get_transcripts(gene_id, gff_db)
    else:
        transcript_ids = list()

    list_annotated_records = list()
    for record in dbnsfp_file.fetch(chrom, start, end):
        record_dict = parse(record, columns_indices, gene_id, transcript_ids)
        if record_dict:
            list_annotated_records.append(record_dict)

    return pd.DataFrame(list_annotated_records)


def parse(record, columns_indices, gene_id, transcript_ids_gff=list()):
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

    return max_value([list_scores[i] for i in list_idx])


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

    return max_value([list_scores[i] for i in list_idx])


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
        max_score = "."

    return max_score


# def max_value(list_scores):
#     """
#     Return the maximum value in a list that can contains "."

#     Args:
#         list_scores (list): list of transcript-based scores for a given metric

#     Returns:
#         float: maximum score in the list
#     """

#     max_score = -100000
#     max_idx = -1
#     for idx, score in enumerate(list_scores):
#         try:
#             if float(score) > max_score:
#                 max_score = score
#                 max_idx = idx
#         except ValueError:
#             continue

#     return max_idx, max_score


def get_transcripts(gene_id, gff_db):
    """
    Return transcripts identifiers of the transcripts found in the GFF file for a given gene

    Args:
        gene_id (str): gene identifier
        gff_db (gffutils.FeatureDB): gffutils database

    Returns:
        list: list of transcripts identifiers
    """
    gene = gff_db[gene_id]

    list_transcript_ids = list()
    for transcript in gff_db.children(gene, featuretype="transcript", order_by="start"):
        list_transcript_ids += [x.split(".")[0] for x in transcript["transcript_id"]]

    return list_transcript_ids


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

    return list_indices, found_column


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


def load_gff(gff_file):
    """Create the gff database used by gffutils.

    Args:
        gff_file (str): Path to a GFF or a gffutils database file.
    Returns:
        gffutils.db: GFF database.

    """

    # gffutils db input
    if gff_file.endswith(".db"):
        gff_db = gffutils.FeatureDB(gff_file)
    # GFF input
    else:
        gff_db_path = gff_file + ".db"
        try:
            os.remove(gff_db_path)
        except OSError:
            pass

        gff_db = gffutils.create_db(gff_file, gff_db_path, merge_strategy="create_unique")

    return gff_db


def remove_duplicates(block_dbnsfp):
    """_summary_

    Args:
        block_dbnsfp (_type_): _description_
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


@click.command()
@click.argument("rates_dnm")
@click.argument("dbnsfp")
@click.argument("annotation_names")
@click.argument("output")
@click.option("--gff", default="")
def annotate_dbnsfp(rates_dnm, dbnsfp, annotation_names, output, gff):
    """
    Annotatate a rates/DNM file with CEP scores from dbNSFP

    Args:
        df (str): Path to rates/DNM file
        dbnsfp (str): Path to dbNSFP file
        annotation_names (str): Path to a file listing which annotations to extract from dbNSFP
        output (str): Path to output file (merged dataframe)
        gff(str) : Path to gff file or gffutils database
    """

    # Load rates/DNM file
    df, add_chr = load_rates_file(rates_dnm, output)
    rates_df_columns = list(df.columns)

    # Load dbnsfp file
    dbnsfp_df = pysam.TabixFile(dbnsfp, encoding="utf-8")

    # Load gff file/DB (to match annotation based on user selected transcripts)
    if gff:
        gff_db = load_gff(gff)
    else:
        gff_db = None

    # Retrieve columns to extract from dbnsfp_df
    dbnsfp_columns_indices, dbnsfp_columns_names = extract_columns_from_dbNFP(dbnsfp_df, annotation_names)

    # For each gene
    list_merged_df = list()
    for gene_id, gene_rates_df in df.groupby("gene_id"):

        chrom = str(gene_rates_df.chrom.values[0]).replace("chr", "")

        # Split each gene in contiguous block (i.e. exons) and load dbNSFP scores
        list_block_df = list()
        for _, block in groupby(sorted(set(gene_rates_df["pos"])), key=lambda n, c=count(): n - next(c)):
            start, end = as_range(list(block))
            try:
                block_dbnsfp = load_dbnsfp(dbnsfp_df, chrom, start - 1, end, dbnsfp_columns_indices, gene_id, gff_db)
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
