import logging.config
import gffutils
import os


def init_log():
    """Initialize logging configuration"""

    MY_LOGGING_CONFIG = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "default_formatter": {"format": "%(levelname)s | %(message)s"},
        },
        "handlers": {
            "stream_handler": {
                "class": "logging.StreamHandler",
                "formatter": "default_formatter",
            },
        },
        "loggers": {
            "logger": {
                "handlers": ["stream_handler"],
                "level": "INFO",
                "propagate": True,
            }
        },
    }

    logging.config.dictConfig(MY_LOGGING_CONFIG)


def read_columns_from_file(columns_file):
    """
    Read columns to extract from annotation file

    Args:
        columns_file (str): file listing columns to extract from TSV/VCF/dbNSFP
    """

    with open(columns_file, "r") as f:
        columns = [x.strip() for x in f.readlines()]

    return columns


def load_gff(gff_file):
    """
    Create the gff database used by gffutils.

    Args:
        gff_file (str): Path to a GFF or a gffutils database file.
    Returns:
        gffutils.db: GFF database.

    """

    logger = logging.getLogger("logger")

    # gffutils db input
    if gff_file.endswith(".db"):
        logger.info(f"Loading GFF db {gff_file}")
        gff_db = gffutils.FeatureDB(gff_file)
    # GFF input
    else:
        gff_db_path = gff_file + ".db"
        logger.info(f"Creating GFF db {gff_db_path}")

        try:
            os.remove(gff_db_path)
        except OSError:
            pass

        gff_db = gffutils.create_db(gff_file, gff_db_path, merge_strategy="create_unique")

    return gff_db
