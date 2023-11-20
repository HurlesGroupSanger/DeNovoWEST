import logging.config
import gffutils
import os

CDS_OFFSET = 50


def init_log():
    """Initialize logging configuration"""

    MY_LOGGING_CONFIG = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "default_formatter": {"format": "[%(levelname)s:%(asctime)s] %(message)s"},
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


def load_gff(gff_file):
    """
    Load a GFF file or a gffutils database file.

    Args:
        gff_file (str): Path to a GFF or a gffutils database file.

    Returns:
        gffutils.db: GFF database.

    """

    logger = logging.getLogger("logger")

    # If gff_file points already to a gffutils database we just load it
    if gff_file.endswith(".db"):
        logger.info(f"Loading gffutils database : {gff_file}")
        gff_db = gffutils.FeatureDB(gff_file)
    # If it points to a GFF file we create the gffutils database
    else:
        gff_db_path = gff_file.replace(".gz", "").replace(".gff", ".gffutils.db")
        logger.info(f"Creating gffutil database : {gff_db_path}")
        gff_db = gffutils.create_db(gff_file, gff_db_path, merge_strategy="create_unique")

    return gff_db
