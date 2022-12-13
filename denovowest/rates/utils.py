import logging.config
import sys
from pathlib import Path

import yaml


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
    """ Superseed configuration in config file with parameters from the command line

    Args:
        conf (dict): configuration from the config file (if any)
        command_params (dict): command line parameters

    Returns:
        dict: Configuration dictionnary.
    """

    for key, value in command_params.items() :
        if value and key != "config":
            conf[key.upper()] = value

    return conf

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
