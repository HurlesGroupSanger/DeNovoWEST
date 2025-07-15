import logging.config


def init_log(show_time=True):
    """
    Initialise logger

    Args:
        show_time (bool, optional): print time in log entries. Defaults to True.
    """

    if show_time:
        log_format = "[%(levelname)s:%(asctime)s] %(message)s"
    else:
        log_format = "%(levelname)s | %(message)s"

    MY_LOGGING_CONFIG = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "default_formatter": {"format": log_format},
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
