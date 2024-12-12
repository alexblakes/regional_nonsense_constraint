""""Package docstring."""

import logging
import os
from pathlib import Path
import sys


def get_module_dot_path():
    """Get the current module's path in dot notation."""
    cwd = os.getcwd()
    abs_path = Path(sys.argv[0])
    rel_path = abs_path.relative_to(cwd)
    rel_path_no_suffix = rel_path.with_suffix("")
    return ".".join(rel_path_no_suffix.parts)


def get_log_file_path(dir="data/logs"):
    if not os.path.exists(dir):
        raise FileNotFoundError(f'Directory "{dir}" does not exist.')

    cwd = os.getcwd()
    rel_path = Path(sys.argv[0]).relative_to(cwd)
    dot_stem = str(rel_path.parent).replace(os.sep, ".")
    log_name = rel_path.with_suffix(".log").name
    log_file = ".".join([dot_stem, log_name])

    return Path(dir).joinpath(log_file)


def setup_logger(name="src", level=logging.DEBUG, stream=True, log_file=True):

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False  # Disconnect from the root logger

    if logger.hasHandlers():
        return logger  # Imported modules won't create a separate log file

    formatter = logging.Formatter(
        "{asctime}|{levelname}|{name}.{funcName}|#{lineno:d}\n{message}",
        style="{",
        datefmt="%d-%m-%y %H:%M:%S",
    )

    if stream:
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(level)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    if log_file:
        if type(log_file) is bool:
            log_file_path = get_log_file_path()  # The default path to log file
        else:
            log_file_path = Path(log_file)  # Allows a bespoke log file

        logger.info(
            f"Writing log to: {log_file_path}"
        )  # Stream handler already defined

        file_handler = logging.FileHandler(log_file_path, mode="w")
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


# A sponge for loggers from external modules
root_logger = setup_logger(name="", log_file=None, stream=True)
