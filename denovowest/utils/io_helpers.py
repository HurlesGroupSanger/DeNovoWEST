"""Shared I/O helpers used by multiple DeNovoWEST commands.

These helpers intentionally stay lightweight: they mostly wrap common file-loading
patterns so command modules can share consistent handling of configuration files,
annotation column lists, and GFF databases.
"""

import logging
from pathlib import Path

import gffutils


def read_columns_from_file(columns_file):
    """Read one annotation column name per line from a plain-text file.

    Args:
        columns_file (str): File listing columns to extract from TSV, VCF, or dbNSFP resources.
    """

    with open(columns_file, "r") as f:
        columns = [x.strip() for x in f.readlines()]

    return columns


def load_gff(gff_file: Path, gff_db_out: Path = None):
    """Create or load a ``gffutils`` database.

    The helper accepts either an existing ``.db`` file or a source GFF file. For
    GFF inputs, the current behaviour is to rebuild the database at the requested
    output path, replacing any existing file at that location.

    Args:
        gff_file (Path): Path to a GFF file or an existing ``gffutils`` database.
        gff_db_out (Path): Optional output path for the created database when ``gff_file`` is a GFF.

    Returns:
        gffutils.FeatureDB: Loaded or newly created database.
    """

    logger = logging.getLogger("logger")
    gff_file = Path(gff_file)

    # gffutils db input
    if gff_file.suffix == ".db":
        logger.info(f"Loading gffutils database {gff_file}")
        gff_db = gffutils.FeatureDB(gff_file)
    # GFF input
    else:
        gff_db_path = Path(gff_db_out) if gff_db_out else gff_file.with_suffix(".db")
        logger.info(f"Creating GFF db {gff_db_path}")

        try:
            gff_db_path.unlink()
            logger.info(f"Removed old gffutils database : {gff_db_path}")
        except OSError:
            pass

        logger.info(f"Creating gffutils database : {gff_db_path}")
        gff_db = gffutils.create_db(str(gff_file), str(gff_db_path), merge_strategy="create_unique")

    return gff_db
