#!/usr/bin/env python

import logging.config


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

consequences_severities = {
    "intergenic": 1,
    "feature_truncation": 2,
    "feature_elongation": 2,
    "regulatory": 3,
    "TF_binding_site": 4,
    "TFBS": 4,
    "downstream": 5,
    "upstream": 5,
    "non_coding_transcript": 6,
    "non_coding": 6,
    "intron": 7,
    "NMD_transcript": 7,
    "non_coding_transcript_exon": 8,
    "5_prime_utr": 9,
    "3_prime_utr": 9,
    "coding_sequence": 10,
    "mature_miRNA": 10,
    "stop_retained": 11,
    "start_retained": 11,
    "synonymous": 11,
    "incomplete_terminal_codon": 12,
    "splice_region": 13,
    "missense": 14,
    "inframe": 14,
    "inframe_insertion": 14,
    "inframe_deletion": 14,
    "protein_altering": 14,
    "transcript_amplification": 15,
    "exon_loss": 16,
    "disruptive": 17,
    "start_lost": 18,
    "stop_lost": 18,
    "stop_gained": 18,
    "frameshift": 18,
    "splice_acceptor": 19,
    "splice_donor": 19,
    "transcript_ablation": 20,
}


DEFAULT_MAX_NB_MUTATIONS_SIM = 250
DEFAULT_MIN_NB_SIM = 10**7

# These inframe/missense and frameshift/nonsense ratios were observed on several databases
# such ad EXAC, gnomAD, TOPMED or DDD.
INFRAME_MISSENSE_RATIO = 0.03
FRAMESHIFT_NONSENSE_RATIO = 1.3

JOBS = 1
RUNTYPE = "ns"


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

    global JOBS
    JOBS = conf["jobs"]

    global RUNTYPE
    RUNTYPE = conf["runtype"]

    logger = logging.getLogger("logger")
    logger.info("DNW simulation running with the following parameters :")
    logger.info("----------")
    for key, value in sorted(conf.items()):
        logger.info(f"{key} : {value}")
    logger.info("----------")
