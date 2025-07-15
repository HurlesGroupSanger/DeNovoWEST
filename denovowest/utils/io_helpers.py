import os
import gffutils
import logging
import yaml


def read_columns_from_file(columns_file):
    """
    Read columns to extract from annotation file

    Args:
        columns_file (str): file listing columns to extract from TSV/VCF/dbNSFP
    """

    with open(columns_file, "r") as f:
        columns = [x.strip() for x in f.readlines()]

    return columns


def load_gff(gff_file, gff_db_out=""):
    """
    Create or load a gffutils database.

    Args:
        gff_file (str): Path to a GFF or a gffutils database file.
        gff_db_out (str): If a GFF file is provided, the user can name the yet to be created gffutils db file.
    Returns:
        gffutils.db: GFF database.

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
    """Load the content of a YAML configuration file in a dictionnary

    Args:
        filename (str): Path to a YAML configuration file.

    Returns:
        dict: Configuration dictionnary.
    """

    with open(filename) as file:
        conf = yaml.load(file, Loader=yaml.FullLoader)

    return conf


def superseed_conf(conf, command_params):
    """Superseed configuration in config file with parameters from the command line

    Args:
        conf (dict): configuration from the config file (if any)
        command_params (dict): command line parameters

    Returns:
        dict: Configuration dictionnary.
    """

    for key, value in command_params.items():
        if value and key != "config":
            conf[key.upper()] = value

    return conf
