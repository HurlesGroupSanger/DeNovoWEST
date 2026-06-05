"""Shared I/O helpers used by multiple DeNovoWEST commands.

These helpers intentionally stay lightweight: they mostly wrap common file-loading
patterns so command modules can share consistent handling of configuration files,
annotation column lists, and GFF databases.
"""

import logging
from pathlib import Path
from typing import Optional

import gffutils
import pandas as pd
import yaml


###################
# File utilities  #
###################


def read_columns_from_file(columns_file: str) -> list[str]:
    """Read one annotation column name per line from a plain-text file.

    Args:
        columns_file: File listing columns to extract from TSV, VCF, or dbNSFP resources.

    Returns:
        List of column names.
    """

    with open(columns_file, "r") as f:
        return [line.strip() for line in f]


def load_conf(filename: str) -> dict:
    """Load the content of a YAML configuration file into a dictionary.

    Args:
        filename: Path to a YAML configuration file.

    Returns:
        Configuration dictionary.
    """

    with open(filename) as f:
        return yaml.load(f, Loader=yaml.FullLoader)


def superseed_conf(conf: dict, command_params: dict) -> dict:
    """Override configuration file values with command-line parameters.

    Args:
        conf: Configuration from the config file.
        command_params: Command-line parameters.

    Returns:
        Updated configuration dictionary.
    """

    for key, value in command_params.items():
        if value and key != "config":
            conf[key.upper()] = value
    return conf


###################
# GFF utilities   #
###################


def load_gff(gff_file: str, gff_db_out: Optional[str] = None) -> gffutils.FeatureDB:
    """Create or load a gffutils database.

    The helper accepts either an existing ``.db`` file or a source GFF file. For
    GFF inputs, the database is rebuilt at the requested output path, replacing
    any existing file at that location.

    Args:
        gff_file: Path to a GFF file or an existing gffutils database.
        gff_db_out: Optional output path for the created database when ``gff_file`` is a GFF.

    Returns:
        Loaded or newly created database.
    """

    logger = logging.getLogger("logger")
    path = Path(gff_file)

    if path.suffix == ".db":
        logger.info(f"Loading gffutils database {gff_file}")
        return gffutils.FeatureDB(str(path))

    db_path = Path(gff_db_out) if gff_db_out else path.with_suffix(".db")
    logger.info(f"Creating GFF db {db_path}")

    if db_path.exists():
        db_path.unlink()
        logger.info(f"Removed old gffutils database: {db_path}")

    return gffutils.create_db(str(path), str(db_path), merge_strategy="create_unique")


def extract_ensembl_gene_id_without_version(gff_db: gffutils.FeatureDB) -> dict[str, str]:
    """Build a mapping from versionless ENSEMBL gene ID to the versioned form stored in the GFF.

    Args:
        gff_db: gffutils database.

    Returns:
        Dictionary mapping bare ENSG ID to versioned ENSG ID
        (e.g. ``{"ENSG00000010404": "ENSG00000010404.1"}``).
    """

    return {
        gene.id.split(".")[0]: gene.id
        for gene in gff_db.features_of_type("gene", order_by="start")
    }


###################
# Genomic helpers #
###################


def as_range(region) -> tuple[int, int]:
    """Return the first and last positions of a genomic block.

    Args:
        region: Iterable of genomic positions.

    Returns:
        Tuple of (start, end) positions.
    """

    positions = list(region)
    return positions[0], positions[-1]


def is_chr_prefixed(file_path: str) -> bool:
    """Return True if chromosome identifiers in the file are prefixed with "chr".

    Inspects up to 100 data rows, looking for a recognised chromosome column
    (``chrom``, ``chr``, ``#chr``, ``#chrom``). Falls back to the first column
    if none of those names are present.

    Args:
        file_path: Path to a tab-separated variant or annotation file.

    Returns:
        True if chromosomes use the "chr" prefix, False otherwise.
    """

    df = pd.read_csv(file_path, sep="\t", comment="#", nrows=100)
    for key in ["chrom", "chr", "#chr", "#chrom"]:
        if key in df.columns:
            return str(df[key].iloc[0]).startswith("chr")
    return str(df.iloc[0, 0]).startswith("chr")
