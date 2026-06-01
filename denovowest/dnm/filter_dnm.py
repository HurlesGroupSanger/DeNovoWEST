#!/usr/bin/env python
import logging
import re
from pathlib import Path

import click
import gffutils
import numpy as np
import pandas as pd

from denovowest.utils.io_helpers import load_gff
from denovowest.utils.log import init_log
from denovowest.utils.params import CDS_OFFSET

#############################
# CLI                       #
#############################


@click.command()
@click.argument("dnm", type=click.Path(exists=True, path_type=Path))
@click.argument("gene_list", type=click.Path(exists=True, path_type=Path))
@click.option("--output_kept_dnm", type=click.Path(path_type=Path))
@click.option("--output_discarded_dnm", type=click.Path(path_type=Path))
@click.option("--gff", type=click.Path(exists=True, path_type=Path))
def filter_dnm(
    dnm: Path,
    gene_list: Path,
    output_kept_dnm: Path,
    output_discarded_dnm: Path,
    gff: Path | None,
) -> None:
    """Remove from the DNM table variants whose gene is absent from the gene list.

    The gene list is either user-provided or built from the rates file.
    If ``--gff`` is provided, variants outside CDS regions are also discarded.

    Args:
        dnm: Path to the DNM file.
        gene_list: Path to the gene list file.
        output_kept_dnm: Destination path for retained variants.
        output_discarded_dnm: Destination path for discarded variants.
        gff: Path to a GFF file or gffutils database (enables CDS filtering).
    """

    init_log()

    genes = load_gene_list(gene_list)
    dnm_df = load_dnm(dnm, genes)

    # Filter on gene list membership
    dnm_df, dnm_discarded_df = filter_on_gene_list(dnm_df, genes)

    # Optionally filter on CDS regions
    if gff:
        dnm_df, dnm_discarded_gff_df = filter_on_gff(dnm_df, gff)
        dnm_discarded_df = pd.concat([dnm_discarded_df, dnm_discarded_gff_df])

    export_dnm(dnm_df, dnm_discarded_df, output_kept_dnm, output_discarded_dnm)
    log_stats(dnm_df, dnm_discarded_df, output_discarded_dnm)


#############################
# Filtering                 #
#############################


def filter_on_gene_list(
    dnm_df: pd.DataFrame,
    gene_list: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Remove rows whose gene is not in ``gene_list``.

    Args:
        dnm_df: DNM table.
        gene_list: Allowed gene identifiers.

    Returns:
        Tuple of (kept DataFrame, discarded DataFrame).
        The discarded DataFrame carries a ``reason`` column.
    """

    kept = dnm_df.loc[dnm_df["gene_id"].isin(gene_list)].copy()
    discarded = dnm_df.loc[~dnm_df["gene_id"].isin(gene_list)].copy()
    discarded["reason"] = "not_in_gene_list"
    return kept, discarded


def filter_on_gff(
    dnm_df: pd.DataFrame,
    gff: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Remove rows that fall outside CDS regions defined in the GFF file.

    Args:
        dnm_df: DNM table (already filtered by gene list).
        gff: Path to a GFF file or gffutils database.

    Returns:
        Tuple of (kept DataFrame, discarded DataFrame).
        The discarded DataFrame carries a ``reason`` column.
    """

    logger = logging.getLogger("logger")
    gff_db = load_gff(gff)

    gene_ids = list(dnm_df["gene_id"].unique())
    genes_cds = build_gene_cds_intervals(gff_db, gene_ids, cds_offset=CDS_OFFSET)
    logger.info("Built gene CDS intervals from GFF")

    list_new_dnm: list[pd.DataFrame] = []
    for gene_id, gene_dnm_df in dnm_df.groupby("gene_id"):
        gene_dnm_df = gene_dnm_df.copy()
        gene_cds = genes_cds[gene_id]

        # Gene absent from the GFF — discard all its variants
        if gene_cds is None:
            logger.warning(f"{gene_id} not found in GFF")
            gene_dnm_df["in_cds"] = False
            gene_dnm_df["reason"] = "gene_not_in_gff"
            list_new_dnm.append(gene_dnm_df)
            continue

        # Vectorised CDS membership check
        hit = np.zeros(len(gene_dnm_df), dtype=bool)
        for cds_start, cds_end in gene_cds:
            hit |= (gene_dnm_df["pos"].values >= cds_start) & (gene_dnm_df["pos"].values <= cds_end)

        gene_dnm_df["in_cds"] = hit
        gene_dnm_df["reason"] = gene_dnm_df["in_cds"].apply(lambda x: "not_in_cds" if not x else "")
        list_new_dnm.append(gene_dnm_df)

    new_dnm_df = pd.concat(list_new_dnm)
    kept = new_dnm_df.loc[new_dnm_df["in_cds"]].copy().drop(columns=["in_cds", "reason"])
    discarded = new_dnm_df.loc[~new_dnm_df["in_cds"]].copy().drop(columns=["in_cds"])
    return kept, discarded


def build_gene_cds_intervals(
    gff_db: gffutils.FeatureDB,
    gene_ids: list[str],
    cds_offset: int,
) -> dict[str, list[tuple[int, int]] | None]:
    """Build CDS intervals for each gene, extended by ``cds_offset`` on each side.

    Args:
        gff_db: GFF utils database.
        gene_ids: Gene identifiers to query.
        cds_offset: Number of extra bases added on each side of every CDS.

    Returns:
        Mapping from gene ID to a list of ``(start, end)`` intervals.
        Genes absent from the GFF map to ``None``.
    """

    intervals: dict[str, list[tuple[int, int]] | None] = {}
    for gid in gene_ids:
        try:
            gene = gff_db[gid]
        except gffutils.exceptions.FeatureNotFoundError:
            intervals[gid] = None
            continue

        raw: list[tuple[int, int]] = []
        for transcript in gff_db.children(gene, level=1):
            for cds in gff_db.children(transcript, featuretype="CDS", order_by="start"):
                start = max(1, cds.start - cds_offset)
                end = cds.end + cds_offset
                raw.append((start, end))

        intervals[gid] = raw if raw else []

    return intervals


#############################
# Data loading / formatting #
#############################


def load_dnm(dnm: Path, genes: list[str]) -> pd.DataFrame:
    """Load the DNM file and normalise gene identifiers.

    Args:
        dnm: Path to the DNM file.
        genes: Gene identifiers from the rates/gene-list file,
            used to reconcile ENSEMBL version differences.

    Returns:
        DNM DataFrame with normalised ``gene_id`` column.
    """

    dnm_df = pd.read_csv(dnm, sep="\t")
    dnm_df["gene_id"] = format_gene_id(list(dnm_df["gene_id"]), genes)
    return dnm_df


def load_gene_list(gene_list: Path) -> list[str]:
    """Load the gene list used to filter the DNM table.

    Args:
        gene_list: Path to the gene list file (one identifier per line).

    Returns:
        List of stripped gene identifiers.
    """

    with open(gene_list, "r") as f:
        return [line.strip() for line in f if line.strip()]


def format_gene_id(genes_in_dnm: list[str], genes_list: list[str]) -> list[str]:
    """Normalise gene identifiers so DNM and rates files use the same format.

    Handles the common mismatch where the rates file has versioned ENSEMBL IDs
    (e.g. ``ENSG00000012048.19``) while the DNM file does not (e.g. ``ENSG00000012048``).

    Args:
        genes_in_dnm: Gene identifiers from the DNM file.
        genes_list: Gene identifiers from the rates/gene-list file.

    Returns:
        Normalised gene identifier list for the DNM file.
    """

    # Strip leading/trailing whitespace
    genes_in_dnm = [g.strip() for g in genes_in_dnm]

    if not genes_in_dnm or not genes_list:
        return genes_in_dnm

    pattern_versioned = r"^ENSG\d{11}\.\d+$"
    pattern_unversioned = r"^ENSG\d{11}$"

    dnm_unversioned = bool(re.match(pattern_unversioned, genes_in_dnm[0]))
    list_versioned = bool(re.match(pattern_versioned, genes_list[0]))

    # If DNM file lacks version but rates file has it, append the version
    if dnm_unversioned and list_versioned:
        version_map = {x.split(".")[0]: x for x in genes_list}
        genes_in_dnm = [version_map.get(g, g) for g in genes_in_dnm]

    return genes_in_dnm


#############################
# Export / reporting        #
#############################


def export_dnm(
    dnm_df: pd.DataFrame,
    dnm_discarded_df: pd.DataFrame,
    output_kept_dnm: Path,
    output_discarded_dnm: Path,
) -> None:
    """Write kept and discarded DNM tables to TSV files.

    Args:
        dnm_df: Retained DNM rows.
        dnm_discarded_df: Discarded DNM rows.
        output_kept_dnm: Destination path for retained variants.
        output_discarded_dnm: Destination path for discarded variants.
    """

    dnm_df.sort_values(by=["chrom", "pos", "ref", "alt"], inplace=True)
    dnm_df.to_csv(output_kept_dnm, sep="\t", index=False)

    dnm_discarded_df.sort_values(by=["chrom", "pos", "ref", "alt"], inplace=True)
    dnm_discarded_df.to_csv(output_discarded_dnm, sep="\t", index=False)


def log_stats(
    dnm_kept_df: pd.DataFrame,
    dnm_discarded_df: pd.DataFrame,
    output_discarded_dnm: Path,
) -> None:
    """Log filtering statistics.

    Args:
        dnm_kept_df: Retained DNM rows.
        dnm_discarded_df: Discarded DNM rows.
        output_discarded_dnm: Path to the discarded DNM file (for user reference).
    """

    logger = logging.getLogger("logger")
    nb_dnm = len(dnm_kept_df) + len(dnm_discarded_df)

    if dnm_discarded_df.empty:
        logger.info(f"All DNM ({nb_dnm}) were retained")
    else:
        nb_discarded = len(dnm_discarded_df)
        logger.warning(f"{nb_discarded}/{nb_dnm} DNM were discarded, of which :")
        for reason, count in dnm_discarded_df["reason"].value_counts().items():
            logger.warning(f"- {reason} : {count} DNM")
        logger.warning(f"Check filtered DNM table : {output_discarded_dnm}")


if __name__ == "__main__":
    filter_dnm()
