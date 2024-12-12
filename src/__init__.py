""""Package docstring."""

import logging
import os
from pathlib import Path
import sys

LOGGER_NAME = "src"
LOG_FORMAT = logging.Formatter(
    "{asctime}|{levelname}|{module}.{funcName}|#{lineno:d} - {message}",
    style="{",
    datefmt="%d-%m-%y %H:%M:%S",
)


def get_module_dot_path():
    """Get the current module's path in dot notation."""
    cwd = os.getcwd()
    abs_path = Path(sys.argv[0])
    rel_path = abs_path.relative_to(cwd)
    rel_path_no_suffix = rel_path.with_suffix("")
    return ".".join(rel_path_no_suffix.parts)


def add_stream_handler(logger, level):
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(LOG_FORMAT)
    logger.addHandler(stream_handler)
    pass


def get_log_file_path(dir="data/logs"):
    if not os.path.exists(dir):
        raise FileNotFoundError(f'Directory "{dir}" does not exist.')

    cwd = os.getcwd()
    rel_path = Path(sys.argv[0]).relative_to(cwd)
    dot_stem = str(rel_path.parent).replace(os.sep, ".")
    log_name = rel_path.with_suffix(".log").name
    log_file = ".".join([dot_stem, log_name])

    return Path(dir).joinpath(log_file)


def setup_logger(name=LOGGER_NAME, level=logging.DEBUG, stream=True):
    """Create a package-level logger.

    Only a stream handler is included by default. Logs are not written to file unless a
    file handler is added separately.
    """

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False  # Disconnect from the root logger

    if logger.hasHandlers():
        return logger
    
    if stream:
        add_stream_handler(logger, level)

    return logger


def add_log_handlers(name=LOGGER_NAME, level=logging.DEBUG, stream=True, log_file=True):
    """Add handlers to the package-level logger.

    Imported modules may call the package logger before the __main__ module does. We
    want only handlers from the __main__ module to persist. This function should be
    called from the __main__ module, after the if __name__ == "__main__" clause.
    """

    logger = logging.getLogger(name)

    logger.handlers.clear()

    if stream:
        add_stream_handler(logger, level)

    if not log_file:
        return logger
    elif log_file is True:  # Use the default log file
        log_file_path = get_log_file_path()
    elif isinstance(log_file, str):  # Specify a bespoke log file
        log_file_path = log_file  # Valid path checks are done behind the scenes
    else:
        raise ValueError(
            f'Invalid parameter for log_file: "{log_file}". '
            f"Must be of type bool or str."
        )

    logger.info(f"Logging to file: {log_file_path}")  # Stream handler already in place

    file_handler = logging.FileHandler(log_file_path, mode="w")
    file_handler.setLevel(level)
    file_handler.setFormatter(LOG_FORMAT)
    logger.addHandler(file_handler)

    return logger


# The root logger acts as a sponge for loggers from external modules:
root_logger = setup_logger(name="", stream=True)

# This is the main logger for the package:
logger = setup_logger()
