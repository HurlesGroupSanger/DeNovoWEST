#!/usr/bin/env python

import logging.config


CONSEQUENCES_MAPPING = {
    "frameshift_variant": "frameshift",
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
    "conserved_exon_terminus_variant": "splice_lof",
}

N_MUTATION_MAX_SIM = 250


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


def log_configuration(conf):
    """_summary_

    Args:
        conf (_type_): _description_
    """

    logger = logging.getLogger("logger")
    logger.info("DNW simulation running with the following parameters :")
    logger.info("----------")
    for key, value in sorted(conf.items()):
        logger.info(f"{key} : {value}")
    logger.info("----------")
