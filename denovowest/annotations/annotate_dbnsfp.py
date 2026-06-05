#!/usr/bin/env python
import logging
import sys
from itertools import count, groupby
from typing import Optional

import click
import pandas as pd
import pysam

from denovowest.utils.io_helpers import (
    as_range,
    extract_ensembl_gene_id_without_version,
    load_gff,
    read_columns_from_file,
)
from denovowest.utils.log import init_log

# TODO : See if we could improve the selection when multiple records exist for a same variant (@remove_duplicates)


def _get_key_column_indices(dbnsfp: pysam.TabixFile) -> tuple[int, int]:
    """Return the column indices of Ensembl_geneid and Ensembl_transcriptid in a dbNSFP file.

    These columns are used to match records to genes and transcripts and their position
    has changed across dbNSFP versions, so they must be resolved from the header at runtime.

    Args:
        dbnsfp: dbNSFP tabix file.

    Returns:
        Tuple of (ensembl_geneid_idx, ensembl_transcriptid_idx).

    Raises:
        ValueError: If either column is absent from the header.
    """

    cols = dbnsfp.header[0].split("\t")
    try:
        return cols.index("Ensembl_geneid"), cols.index("Ensembl_transcriptid")
    except ValueError as e:
        raise ValueError(f"Required dbNSFP column not found in header: {e}") from e


def load_dbnsfp(
    dbnsfp_file: pysam.TabixFile,
    chrom: str,
    start: int,
    end: int,
    columns_indices: list[int],
    gene_id: str,
    gff_db=None,
    ensembl_gene_id_map_version: Optional[dict[str, str]] = None,
    ensembl_geneid_idx: int = 15,
    ensembl_transcriptid_idx: int = 16,
) -> pd.DataFrame:
    """Fetch dbNSFP records for a genomic region.

    Args:
        dbnsfp_file: dbNSFP tabix file.
        chrom: Chromosome identifier (without "chr" prefix).
        start: 0-based start position.
        end: End position.
        columns_indices: Indices of dbNSFP columns to retrieve.
        gene_id: Gene identifier used in the DNM/rates file.
        gff_db: gffutils database (optional; used to retrieve transcripts).
        ensembl_gene_id_map_version: Maps versionless ENSG to versioned form.
        ensembl_geneid_idx: Column index of Ensembl_geneid (resolved from header).
        ensembl_transcriptid_idx: Column index of Ensembl_transcriptid (resolved from header).

    Returns:
        DataFrame of annotated records.
    """

    # Retrieve the transcripts associated to the gene from the GFF file
    if gff_db:
        transcript_ids = get_transcripts(gene_id, gff_db, ensembl_gene_id_map_version or {})
    else:
        transcript_ids = []

    list_annotated_records = []
    for record in dbnsfp_file.fetch(chrom, start, end):
        record_dict = parse(
            record, columns_indices, gene_id, transcript_ids,
            ensembl_geneid_idx, ensembl_transcriptid_idx,
        )
        if record_dict:
            list_annotated_records.append(record_dict)

    return pd.DataFrame(list_annotated_records)


def get_transcripts(
    gene_id: str,
    gff_db,
    ensembl_gene_id_map_version: dict[str, str],
) -> list[str]:
    """Return transcript identifiers from the GFF file for a given gene.

    Args:
        gene_id: Gene identifier.
        gff_db: gffutils database.
        ensembl_gene_id_map_version: Maps versionless ENSG to versioned form.

    Returns:
        List of transcript identifiers (without version suffix).
    """

    try:
        if "." in gene_id:
            gene = gff_db[gene_id]
        else:
            gene = gff_db[ensembl_gene_id_map_version[gene_id]]

        list_transcript_ids = []
        for transcript in gff_db.children(gene, featuretype="transcript", order_by="start"):
            list_transcript_ids += [x.split(".")[0] for x in transcript["transcript_id"]]
    except KeyError:
        list_transcript_ids = []

    return list_transcript_ids


def parse(
    record: str,
    columns_indices: list[int],
    gene_id: str,
    transcript_ids_gff: Optional[list[str]] = None,
    ensembl_geneid_idx: int = 15,
    ensembl_transcriptid_idx: int = 16,
) -> dict:
    """Parse a dbNSFP record and extract the data of interest.

    Args:
        record: Raw tab-separated line from dbNSFP.
        columns_indices: Indices of columns to retrieve.
        gene_id: ENSEMBL gene identifier from the input file.
        transcript_ids_gff: Transcripts for this gene found in the GFF (if provided).
        ensembl_geneid_idx: Column index of Ensembl_geneid (version-dependent; use
            _get_key_column_indices to resolve from the file header).
        ensembl_transcriptid_idx: Column index of Ensembl_transcriptid.

    Returns:
        Dictionary with variant coordinates and scores, or empty dict if the record
        belongs to a different (overlapping) gene.
    """

    if transcript_ids_gff is None:
        transcript_ids_gff = []

    fields = record.split("\t")

    record_dict: dict = {}
    record_dict["chrom"] = fields[0]
    record_dict["pos"] = int(fields[1])
    record_dict["ref"] = fields[2]
    record_dict["alt"] = fields[3]

    ensembl_gene_ids_dbnsfp = fields[ensembl_geneid_idx].split(";")
    transcript_ids_dbnsfp = fields[ensembl_transcriptid_idx].split(";")

    # Check the intersection between the transcripts from the user GFF and dbNSFP
    transcript_in_dbnsfp = bool(set(transcript_ids_gff) & set(transcript_ids_dbnsfp))
    record_dict["transcript_in_dbnsfp"] = transcript_in_dbnsfp

    # If the record is associated to an overlapping gene, skip it
    if gene_id.split(".")[0] not in ensembl_gene_ids_dbnsfp:
        return {}

    # Retrieve the scores
    for column_index in columns_indices:
        record_dict[str(column_index)] = fields[column_index]

    # If a GFF has been provided, retrieve the maximum transcript-based score
    if transcript_ids_gff:
        for column_index in columns_indices:
            list_scores = fields[column_index].split(";")

            # If there is only one score, at this point we assume it is a per-gene score
            if len(list_scores) == 1:
                record_dict[str(column_index)] = list_scores[0]
                continue

            score = "."

            # We first try to retrieve the transcript-based score
            if transcript_in_dbnsfp:
                score = get_transcript_based_score(list_scores, transcript_ids_gff, transcript_ids_dbnsfp)

            # Some MANE transcripts are missing in dbNSFP, therefore if no score associated to a transcript
            # could be retrieved we switch to gene mode
            if score == ".":
                score = get_gene_based_score(list_scores, gene_id.split(".")[0], ensembl_gene_ids_dbnsfp)

            record_dict[str(column_index)] = score

    return record_dict


def get_transcript_based_score(
    list_scores: list[str],
    transcript_ids_gff: list[str],
    transcript_ids_dbnsfp: list[str],
) -> str:
    """Return the maximum score among transcripts found in the GFF.

    Args:
        list_scores: Per-transcript scores for a given metric.
        transcript_ids_gff: Transcript identifiers found in the GFF.
        transcript_ids_dbnsfp: Transcript identifiers in dbNSFP for this variant.

    Returns:
        Maximum score string, or the raw semicolon-joined value if scores cannot be parsed.
    """

    list_idx = [idx for idx, tid in enumerate(transcript_ids_dbnsfp) if tid in transcript_ids_gff]

    # Not all dbNSFP columns contain scores, and they do not necessarily match the number of values in gene/transcript columns
    try:
        max_score = max_value([list_scores[i] for i in list_idx])
    except IndexError:
        # In that case we just return the whole value
        max_score = ";".join(list_scores)

    return max_score


def get_gene_based_score(
    list_scores: list[str],
    gene_id: str,
    gene_ids_dbnsfp: list[str],
) -> str:
    """Return the maximum score across all transcripts for this gene in dbNSFP.

    Args:
        list_scores: Per-transcript scores for a given metric.
        gene_id: Versionless ENSEMBL gene identifier.
        gene_ids_dbnsfp: Gene identifiers in dbNSFP for this variant.

    Returns:
        Maximum score string, or the raw semicolon-joined value if scores cannot be parsed.
    """

    list_idx = [idx for idx, gid in enumerate(gene_ids_dbnsfp) if gid == gene_id]

    # Not all dbNSFP columns contain scores, and they do not necessarily match the number of values in gene/transcript columns
    try:
        max_score = max_value([list_scores[i] for i in list_idx])
    except IndexError:
        # In that case we just return the whole value
        max_score = ";".join(list_scores)

    return max_score


def max_value(list_scores: list[str]):
    """Return the maximum numeric value in a list that may contain missing values (".").

    Args:
        list_scores: List of score strings.

    Returns:
        Maximum score as a float, "." if all values are missing, or the
        semicolon-joined raw value for non-numeric columns.
    """

    try:
        return max(float(x) for x in list_scores if x != ".")
    except ValueError:
        # Not all dbNSFP columns contain scores, and they do not necessarily match the number of values in gene/transcript columns
        if set(list_scores) == {"."}:
            # If there is no score we return "."
            return "."
        # If it was not a score column we return the whole value
        return ";".join(list_scores)


def extract_columns_from_dbnsfp(
    dbnsfp: pysam.TabixFile,
    columns: str,
    columns_file: str,
) -> tuple[list[int], list[str]]:
    """Extract the indices and names of columns to use from a dbNSFP file.

    Args:
        dbnsfp: dbNSFP tabix file.
        columns: Comma-separated column names to extract.
        columns_file: File listing column names to extract (one per line).

    Returns:
        Tuple of (column_indices, column_names).
    """

    logger = logging.getLogger("logger")

    # If the user provided the columns to extract in a separate file we read it
    if columns_file:
        columns_to_extract = read_columns_from_file(columns_file)
    # If provided as a string we split it
    elif columns:
        columns_to_extract = columns.split(",")
    # Default: all columns beyond the standard [chrom, pos, ref, alt]
    else:
        columns_to_extract = dbnsfp.header[0].split("\t")[4:]

    dbnsfp_columns = dbnsfp.header[0].split("\t")

    # Find the index of all columns asked by the user
    list_indices = []
    found_columns = []
    for idx, name in enumerate(dbnsfp_columns):
        if name in columns_to_extract:
            list_indices.append(idx)
            found_columns.append(name)

    columns_not_found = set(columns_to_extract) - set(found_columns)
    if columns_not_found:
        logger.error(f"Columns not found in dbNSFP: {columns_not_found}")
        sys.exit(1)

    return list_indices, found_columns


def check_columns(input_columns: list[str], dbnsfp_columns: list[str]) -> set[str]:
    """Identify dbNSFP columns that already exist in the variants file.

    Args:
        input_columns: Columns already present in the input variants DataFrame.
        dbnsfp_columns: Columns from dbNSFP.

    Returns:
        Set of column names present in both.
    """

    logger = logging.getLogger("logger")

    shared_columns = set(input_columns) & set(dbnsfp_columns)
    for col in shared_columns:
        logger.warning(
            f"{col} already exists in the variant file; "
            "the dbNSFP column will be suffixed with _dbnsfp"
        )
    return shared_columns


def load_rates_file(rates: str) -> tuple[pd.DataFrame, bool]:
    """Load a rates/DNM file.

    Args:
        rates: Path to the variants file in TSV format.

    Returns:
        Tuple of (DataFrame, add_chr) where add_chr indicates whether chromosomes
        in the file use the "chr" prefix.
    """

    logger = logging.getLogger("logger")

    rates_df = pd.read_table(rates, dtype={"chrom": str, "pos": int, "ref": str, "alt": str})

    # Edge case: splitting processes can yield empty input files
    if rates_df.empty:
        logger.warning("Rates file is empty")
        return rates_df, False

    # Depending on the GFF, chromosome can be defined as "chrN" or just "N"
    add_chr = str(rates_df.iloc[0].chrom).startswith("chr")

    return rates_df, add_chr


def remove_duplicates(block_dbnsfp: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate dbNSFP records for the same variant.

    dbNSFP can have two records for the same variant when it leads to different
    amino acids. We prioritise the record that matches a GFF transcript, then the
    one with the fewest missing annotations.

    Args:
        block_dbnsfp: Subset of the annotated DataFrame for a genomic block.

    Returns:
        Deduplicated DataFrame.
    """

    indices_to_remove = []
    duplicated_mask = block_dbnsfp[["chrom", "pos", "ref", "alt"]].duplicated(keep=False)
    for _, min_df in block_dbnsfp.loc[duplicated_mask].groupby(["pos", "alt"]):

        indices = min_df.index
        if min_df["transcript_in_dbnsfp"].sum() != 0:
            min_df = min_df.loc[min_df["transcript_in_dbnsfp"]].copy()

        # Keep the record with the most annotated values
        min_df["nb_missing_values"] = (min_df != ".").sum(axis=1)
        keep_idx = min_df.sort_values("nb_missing_values", ascending=False).index[0]

        indices_to_remove += list(set(indices) - {keep_idx})

    return block_dbnsfp.drop(indices_to_remove)


@click.command()
@click.argument("rates_dnm")
@click.argument("dbnsfp")
@click.argument("output")
@click.option("-c", "--columns", default="", help="Annotations to extract from dbNSFP (comma-separated)")
@click.option("-C", "--columns-file", default="", help="File listing which annotations to extract from dbNSFP")
@click.option("--gff", default="", help="GFF file or gffutils database (for transcript-level score selection)")
def annotate_dbnsfp(rates_dnm: str, dbnsfp: str, output: str, columns: str, columns_file: str, gff: str) -> None:
    """Annotate a rates/DNM file with scores from dbNSFP.

    When multiple scores exist for a variant (e.g. multiple transcripts, overlapping genes),
    the maximum score among transcripts found in the GFF is retrieved. When no GFF is
    provided, the maximum score for the matching gene is used.

    Args:
        rates_dnm: Variants file (e.g. rates, DNM) with columns [gene_id, chrom, pos, ref, alt, ...].
        dbnsfp: dbNSFP genome-wide tabix file.
        output: Output annotated DataFrame.
        columns: Annotations to extract from dbNSFP (comma-separated).
        columns_file: File listing which annotations to extract from dbNSFP.
        gff: GFF file or gffutils database for transcript-level score selection.
    """

    init_log()
    logger = logging.getLogger("logger")

    # Load rates/DNM file
    df, add_chr = load_rates_file(rates_dnm)
    rates_df_columns = list(df.columns)

    # Load dbnsfp file
    dbnsfp_df = pysam.TabixFile(dbnsfp, encoding="utf-8")

    # Load GFF (to match annotation based on user-selected transcripts)
    if gff:
        gff_db = load_gff(gff)
        ensembl_gene_id_map_version = extract_ensembl_gene_id_without_version(gff_db)
    else:
        gff_db = None
        ensembl_gene_id_map_version = {}

    # Resolve gene/transcript key column indices from the header (position varies across dbNSFP versions)
    ensembl_geneid_idx, ensembl_transcriptid_idx = _get_key_column_indices(dbnsfp_df)

    # Retrieve columns to extract from dbNSFP
    dbnsfp_columns_indices, dbnsfp_columns_names = extract_columns_from_dbnsfp(dbnsfp_df, columns, columns_file)

    # Suffix any column names that clash with the input file
    existing_columns = check_columns(rates_df_columns, dbnsfp_columns_names)
    dbnsfp_columns_names = [f"{col}_dbnsfp" if col in existing_columns else col for col in dbnsfp_columns_names]

    # For each gene
    list_merged_df = []
    for gene_id, gene_rates_df in df.groupby("gene_id"):

        logger.info(f"Annotating {gene_id}")

        chrom = str(gene_rates_df.chrom.values[0]).replace("chr", "")

        # Split each gene in contiguous blocks (i.e. exons) and load dbNSFP scores
        list_block_df = []
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
                    ensembl_geneid_idx,
                    ensembl_transcriptid_idx,
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
                gene_dbnsfp_df["chrom"] = "chr" + gene_dbnsfp_df["chrom"]
            merged_gene_df = gene_rates_df.merge(gene_dbnsfp_df, how="left", on=["chrom", "pos", "ref", "alt"])
        else:
            merged_gene_df = gene_rates_df

        list_merged_df.append(merged_gene_df)

    # Combine results for all genes
    if df.empty:
        # When using Roulette we can end up with an empty input file, in that case we just return the empty data frame
        merged_df = df
    else:
        merged_df = pd.concat(list_merged_df)

    # If there is no data in dbNSFP for any gene we fill the columns with NA
    if merged_df.shape[1] == len(rates_df_columns):
        merged_df["transcript_in_dbnsfp"] = False
        for column in dbnsfp_columns_names:
            merged_df[column] = pd.NA
    else:
        merged_df.columns = rates_df_columns + ["transcript_in_dbnsfp"] + dbnsfp_columns_names

    merged_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    annotate_dbnsfp()
