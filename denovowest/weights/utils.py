import logging.config
import logging
import pandas as pd


CONSEQUENCES_MAPPING = {
    "frameshift_variant": "frameshift",
    "frameshift": "frameshift",
    "inframe_insertion": "inframe",
    "inframe_deletion": "inframe",
    "missense_variant": "missense",
    "missense": "missense",
    "stop_gained": "nonsense",
    "synonymous_variant": "synonymous",
    "splice_acceptor_variant": "splice_lof",
    "splice_donor_variant": "splice_lof",
    "splice_acceptor": "splice_lof",
    "splice_donor": "splice_lof",
    "splice_region_variant": "splice_region",
    "splice_region": "splice_region",
    "conserved_exon_terminus_variant": "splice_lof",
    "start_lost": "missense",
    "stop_lost": "missense",
    "stop_retained": "synonymous",
    "synonymous": "synonymous",
    "nonsense": "nonsense",
    "splice_lof": "splice_lof",
    "inframe": "inframe",
}


def filter_on_consequences(df: pd.DataFrame):
    """
    Filter all variants with a consequence not found in CONSEQUENCES_MAPPING

    Args:
        df (pd.DataFrame): variant table (rates or dnm) having a consequence column
    """

    logger = logging.getLogger("logger")

    # TODO : Replace with better consequence extraction
    # When bcftools report several consequences separated by the character "&", we extract the first one
    df.consequence = [csq.split("&")[0] if isinstance(csq, str) else csq for csq in list(df.consequence)]

    filt = df.consequence.isin(CONSEQUENCES_MAPPING.keys())
    kept_df = df.loc[filt]

    logger.info(f"Before consequence filtering : {df.shape[0]} records")

    discarded_df = df.loc[~filt]
    if not discarded_df.empty:
        count_discarded = discarded_df["consequence"].value_counts()
        logger.warning(f"{discarded_df.shape[0]}/{df.shape[0]} records were discarded")
        for consequence, count in count_discarded.items():
            logger.warning(f"{consequence} : {count}")
    else:
        logger.info("All records have an acceptable consequence.")

    logger.info(f"After consequence filtering : {kept_df.shape[0]} records")

    return kept_df


def assign_meta_consequences(df: pd.DataFrame):
    """
    Assign a higher level consequence to each variant

    Args:
        df (pd.DataFrame): variant table (rates or dnm) having a consequence column
    """

    logger = logging.getLogger("logger")

    df.consequence = df.consequence.replace(CONSEQUENCES_MAPPING)

    consequence_counts = dict(df.consequence.value_counts())
    for consequence, count in consequence_counts.items():
        logger.info(f"{consequence} : {count}")

    return df


def init_log():
    """Initialize logging configuration"""

    MY_LOGGING_CONFIG = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "default_formatter": {"format": "[%(levelname)s] %(message)s"},
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
