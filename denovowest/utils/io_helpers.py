"""Shared I/O helpers used by multiple DeNovoWEST commands.

These helpers intentionally stay lightweight: they mostly wrap common file-loading
patterns so command modules can share consistent handling of configuration files,
annotation column lists, and GFF databases.
"""

import os
import gffutils
import logging
import yaml


def read_columns_from_file(columns_file):
    """Read one annotation column name per line from a plain-text file.

    Args:
        columns_file (str): File listing columns to extract from TSV, VCF, or dbNSFP resources.
    """

    with open(columns_file, "r") as f:
        columns = [x.strip() for x in f.readlines()]

    return columns


def load_gff(gff_file, gff_db_out=""):
    """Create or load a ``gffutils`` database.

    The helper accepts either an existing ``.db`` file or a source GFF file. For
    GFF inputs, the current behaviour is to rebuild the database at the requested
    output path, replacing any existing file at that location.

    Args:
        gff_file (str): Path to a GFF file or an existing ``gffutils`` database.
        gff_db_out (str): Optional output path for the created database when ``gff_file`` is a GFF.

    Returns:
        gffutils.FeatureDB: Loaded or newly created database.
    """

    logger = logging.getLogger("logger")

    # gffutils db input
    if gff_file.endswith(".db"):
        logger.info(f"Loading gffutils database {gff_file}")
        gff_db = gffutils.FeatureDB(gff_file)
    # GFF input
    else:
        if gff_db_out:
            gff_db_path = gff_db_out
        else:
            gff_db_path = gff_file + ".db"
        logger.info(f"Creating GFF db {gff_db_path}")

        try:
            os.remove(gff_db_path)
            logger.info(f"Removed old gffutils database : {gff_db_path}")
        except OSError:
            pass

        logger.info(f"Creating gffutils database : {gff_db_path}")
        gff_db = gffutils.create_db(gff_file, gff_db_path, merge_strategy="create_unique")

    return gff_db


def load_conf(filename):
    """Load a YAML configuration file into a Python dictionary.

    Args:
        filename (str): Path to a YAML configuration file.

    Returns:
        dict: Parsed configuration dictionary.
    """

    with open(filename) as file:
        conf = yaml.load(file, Loader=yaml.FullLoader)

    return conf


def superseed_conf(conf, command_params):
    """Override config-file values with non-empty command-line parameters.

    Args:
        conf (dict): Configuration loaded from a file.
        command_params (dict): Parsed command-line parameters.

    Returns:
        dict: Updated configuration dictionary with CLI values taking precedence.
    """

    for key, value in command_params.items():
        if value and key != "config":
            conf[key.upper()] = value

    return conf
