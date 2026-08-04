#!/usr/bin/env python
import gzip
import io
import json
import logging
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Optional

import click
import numpy as np
import pandas as pd
import pysam

from denovowest.utils.log import init_log
from denovowest.utils.params import ConsequenceGroups

###################
# Constants       #
###################

# Structural columns present in every rates file; everything else is an annotation score.
FIXED_COLUMNS: frozenset[str] = frozenset({"gene_id", "chrom", "pos", "ref", "alt", "consequence", "prob"})

FIXED_COLUMN_DTYPES: dict[str, type] = {
    "gene_id": str,
    "chrom": str,
    "pos": int,
    "ref": str,
    "alt": str,
    "consequence": str,
    "prob": float,
}

CHROMOSOMES: list[str] = ["chr" + str(x) for x in range(1, 23)] + ["chrX"]

THRESHOLD_MISSING_ANNOTATIONS: float = 0.1


###################
# Core logic      #
###################


def get_score_columns(df: pd.DataFrame) -> list[str]:
    """Return annotation score column names.

    A column is a score if it is not a fixed structural column AND has float dtype.
    The float check excludes non-score string columns (e.g. variant or transcript IDs).
    Must be called after handle_nas so that score columns have already been coerced.

    Args:
        df: rates DataFrame, after handle_nas has been applied.

    Returns:
        Ordered list of score column names.
    """

    return [col for col in df.columns if col not in FIXED_COLUMNS and pd.api.types.is_float_dtype(df[col])]


def rates_stats(rates_df: pd.DataFrame, deciles: bool = False) -> dict:
    """Compute per-gene annotation statistics.

    Args:
        rates_df: rates DataFrame with at least the FIXED_COLUMNS columns.
        deciles: whether to include score deciles in the output.

    Returns:
        dict mapping gene_id to per-consequence annotation statistics.
    """

    # Compute score columns once to avoid repeating the lookup per group.
    score_cols = get_score_columns(rates_df)
    return {gene_id: _gene_stats(gene_df, score_cols, deciles) for gene_id, gene_df in rates_df.groupby("gene_id")}


def _gene_stats(gene_df: pd.DataFrame, score_cols: list[str], deciles: bool = False) -> dict:
    """Compute per-consequence statistics for a single gene.

    Args:
        gene_df: rates rows for a single gene.
        score_cols: annotation score columns to include in stats.
        deciles: whether to include score deciles in the output.

    Returns:
        dict mapping consequence to annotation statistics.
    """

    return {
        consequence: _consequence_stats(consequence, cons_df, score_cols, deciles)
        for consequence, cons_df in gene_df.groupby("consequence")
    }


def _consequence_stats(consequence: str, consequence_df: pd.DataFrame, score_cols: list[str], deciles: bool = False) -> dict:
    """Compute missing-annotation counts and score distributions for one consequence.

    For missense, score entries are always recorded (including 100% missing) so that
    gaps in missense-only predictors are visible. For all other consequences, entries
    where every variant is NaN are omitted — the predictor simply doesn't apply.

    Args:
        consequence: consequence label (e.g. "missense", "splice_donor").
        consequence_df: rates rows for a single gene/consequence pair.
        score_cols: annotation score columns to include in stats.
        deciles: whether to include score deciles in the output.

    Returns:
        dict with nb_mutations, nb_missing_prob and per-score breakdown.
    """

    nb_mutations = len(consequence_df)
    stats: dict = {
        "nb_mutations": nb_mutations,
        "nb_missing_prob": int(consequence_df["prob"].isna().sum()),
    }

    for score in score_cols:
        nb_missing = int(consequence_df[score].isna().sum())
        # For non-missense consequences, skip scores where every variant is NaN:
        # the predictor doesn't apply (e.g. a missense-only predictor on splice variants).
        if nb_missing == nb_mutations and consequence != "missense":
            continue
        score_stats: dict = {"nb_missing_score": nb_missing}
        if deciles and nb_missing < nb_mutations:
            score_deciles = consequence_df[score].quantile([0.1 * i for i in range(1, 10)])
            score_stats["deciles"] = list(score_deciles)
        stats[score] = score_stats

    return stats


###################
# Helpers         #
###################


def handle_nas(rates_df: pd.DataFrame) -> pd.DataFrame:
    """Replace "." and empty strings with NaN and coerce annotation columns to float.

    Fixed structural columns are skipped; only annotation score columns are coerced.

    Args:
        rates_df: rates DataFrame before NaN normalisation.

    Returns:
        rates DataFrame with NaN-normalised annotation columns.
    """

    rates_df = rates_df.replace({".": np.nan, "": np.nan})
    for col in rates_df.columns:
        if col in FIXED_COLUMNS:
            continue
        try:
            rates_df[col] = rates_df[col].astype(float)
        except (ValueError, TypeError):
            pass
    return rates_df


def warning_genes(rates_stats_dict: dict) -> None:
    """Log a warning for each gene/score pair exceeding the missing-annotation threshold.

    Args:
        rates_stats_dict: output of rates_stats().
    """

    logger = logging.getLogger("logger")
    for gene_id, g_stats in rates_stats_dict.items():
        for score, percentage in percentage_missing_annotations(g_stats).items():
            if percentage > THRESHOLD_MISSING_ANNOTATIONS:
                logger.warning(f"{gene_id} | {score} | {round(percentage * 100, 2)}% missing annotations")


def percentage_missing_annotations(gene_stats_dict: dict) -> dict[str, float]:
    """Compute per-score fraction of variants missing an annotation.

    The denominator for each score only counts consequences where that score has
    an entry: missense-only predictors are therefore evaluated over missense sites
    alone, while all-consequence predictors are evaluated over all scored sites.

    Args:
        gene_stats_dict: per-consequence annotation statistics for one gene.

    Returns:
        dict mapping score name to fraction missing (0–1).
    """

    nb_mutations_per_score: dict[str, int] = {}
    nb_missing: dict[str, int] = {}

    for cons_stats in gene_stats_dict.values():
        nb_cons = cons_stats["nb_mutations"]
        for score, score_stats in cons_stats.items():
            if not isinstance(score_stats, dict):
                continue
            nb_mutations_per_score[score] = nb_mutations_per_score.get(score, 0) + nb_cons
            nb_missing[score] = nb_missing.get(score, 0) + score_stats["nb_missing_score"]

    return {
        score: count / nb_mutations_per_score[score]
        for score, count in nb_missing.items()
        if nb_mutations_per_score.get(score, 0) > 0
    }


###################
# I/O             #
###################


def _read_rates_header(rates_path: str) -> tuple[str, list[str]]:
    """Read the header line from a gzip-compressed rates file.

    Args:
        rates_path: path to the gzip-compressed rates file.

    Returns:
        Tuple of (raw header line, list of column names).
    """

    with gzip.open(rates_path, "rt") as fh:
        header_line = fh.readline().rstrip("\n")
    return header_line, header_line.split("\t")


def _fetch_chromosome_df(
    rates_file: pysam.TabixFile,
    chromosome: str,
    header_line: str,
) -> Optional[pd.DataFrame]:
    """Fetch all rows for one chromosome into a DataFrame via pandas CSV parser.

    Using io.StringIO + pd.read_csv is significantly faster than building a
    DataFrame from a Python-level list of split strings.

    Args:
        rates_file: open pysam TabixFile handle.
        chromosome: chromosome name (e.g. "chr1").
        header_line: raw TSV header line prepended so pandas can parse column names.

    Returns:
        DataFrame for the chromosome, or None if the chromosome is absent or empty.
    """

    try:
        rows = list(rates_file.fetch(chromosome))
    except ValueError:
        return None

    if not rows:
        return None

    content = header_line + "\n" + "\n".join(rows)
    return pd.read_csv(io.StringIO(content), sep="\t", dtype=FIXED_COLUMN_DTYPES, low_memory=False)


def _process_chromosome(rates_path: str, chromosome: str, header_line: str, deciles: bool = False) -> Optional[dict]:
    """Fetch and compute stats for one chromosome.

    Opens its own TabixFile handle so this function is safe to call from threads.

    Args:
        rates_path: path to the gzip-compressed tabix-indexed rates file.
        chromosome: chromosome name (e.g. "chr1").
        header_line: raw TSV header line.
        deciles: whether to include score deciles in the output.

    Returns:
        Per-gene stats dict for the chromosome, or None if absent.
    """

    rates_file = pysam.TabixFile(rates_path, encoding="utf-8")
    chrom_df = _fetch_chromosome_df(rates_file, chromosome, header_line)
    if chrom_df is None:
        return None
    chrom_df = handle_nas(chrom_df)
    chrom_df = chrom_df[chrom_df["consequence"].isin(ConsequenceGroups.OE_RATIO)]
    return rates_stats(chrom_df, deciles=deciles)


###################
# CLI entry point #
###################


@click.command()
@click.argument("rates")
@click.argument("output")
@click.option("--fullgenome", is_flag=True)
@click.option(
    "--threads",
    default=1,
    show_default=True,
    help="Number of chromosomes to process in parallel (--fullgenome only). "
    "Higher values are faster but use proportionally more memory.",
)
@click.option("--deciles", is_flag=True, default=False, help="Include per-score decile distributions in the output (increases file size).")
def main(rates: str, output: str, fullgenome: bool, threads: int, deciles: bool) -> None:
    """Compute per-gene annotation statistics for an annotated mutation rates file.

    Outputs a JSON file with per-gene, per-consequence counts of missing annotations
    and score distributions. Logs warnings for genes with >10% missing annotations.

    Designed to run on subsetted rates files from the Nextflow pipeline. Use
    --fullgenome for whole-genome tabix-indexed files.

    Args:
        rates: mutation rates file (TSV or gzip-compressed tabix-indexed TSV).
        output: JSON file to write per-gene statistics.
        fullgenome: set when rates is a full-genome tabix-indexed file.
        threads: number of chromosomes to process concurrently (fullgenome only).
        deciles: whether to include per-score decile distributions in the output.
    """

    init_log()
    logger = logging.getLogger("logger")

    if fullgenome:
        header_line, _ = _read_rates_header(rates)
        n_workers = min(threads, len(CHROMOSOMES), os.cpu_count() or 1)
        logger.info(f"Processing {len(CHROMOSOMES)} chromosomes with {n_workers} threads")

        # Use threads on a single core (low overhead, GIL released during I/O);
        # use processes on multiple cores to bypass the GIL for CPU-bound pandas work.
        executor_cls = ThreadPoolExecutor if n_workers == 1 else ProcessPoolExecutor

        rates_stats_dict: dict = {}
        with executor_cls(max_workers=n_workers) as pool:
            futures = {pool.submit(_process_chromosome, rates, chrom, header_line, deciles): chrom for chrom in CHROMOSOMES}
            for future in as_completed(futures):
                chrom = futures[future]
                chrom_stats = future.result()
                if chrom_stats is None:
                    logger.warning(f"# {chrom} not found in index")
                else:
                    rates_stats_dict.update(chrom_stats)

    else:
        rates_df = pd.read_csv(rates, sep="\t", dtype=FIXED_COLUMN_DTYPES)
        rates_df = handle_nas(rates_df)
        rates_df = rates_df[rates_df["consequence"].isin(ConsequenceGroups.OE_RATIO)]
        rates_stats_dict = rates_stats(rates_df, deciles=deciles)

    # Export stats in a JSON file
    with open(output, "w") as fh:
        json.dump(rates_stats_dict, fh, indent=2)

    # Log genes with too many missing annotations
    warning_genes(rates_stats_dict)


if __name__ == "__main__":
    main()
