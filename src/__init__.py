""""Package docstring."""

import logging
from pathlib import Path
import sys

from sklego import pandas_utils

LOG_FORMAT = logging.Formatter(
    ">> {asctime}|{levelname}|{module}.{funcName}|#{lineno:d}\n{message}\n",
    style="{",
    datefmt="%d-%m-%y %H:%M:%S",
)


def get_module_dot_path():
    """Get the current module's path in dot notation."""
    cwd = Path.cwd()
    abs_path = Path(sys.argv[0]).resolve()
    rel_path = abs_path.relative_to(cwd)
    rel_path_no_suffix = rel_path.with_suffix("")
    return ".".join(rel_path_no_suffix.parts)


def get_log_file_path(dir="data/logs"):
    dir = Path(dir)
    if not dir.is_dir():
        raise FileNotFoundError(f'Directory "{dir}" does not exist.')

    log_stem = get_module_dot_path()
    log_file = ".".join([log_stem, "log"])

    return dir.joinpath(log_file)


def add_stream_handler(logger, level):
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(LOG_FORMAT)
    logger.addHandler(stream_handler)
    return logger


def add_file_handler(logger, level):
    log_file_path = get_log_file_path()
    file_handler = logging.FileHandler(log_file_path, mode="w")
    file_handler.setLevel(level)
    file_handler.setFormatter(LOG_FORMAT)
    logger.addHandler(file_handler)
    logger.info(f"Logging to file: {log_file_path}")
    return logger


def add_log_handlers(name="src", level=logging.DEBUG, stream=True, file=True):
    logger = logging.getLogger(name)
    logger.handlers.clear()
    if stream:
        add_stream_handler(logger, level)
    if file:
        add_file_handler(logger, level)
    return logger


def setup_logger(name="src", level=logging.DEBUG):
    """Create a package-level logger.

    Only a stream handler is included by default. Logs are not written to file unless a
    file handler is added separately.
    """

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False  # Disconnect from the root logger

    return logger


def write_out(df, path, *args, **kwargs):
    kwargs.setdefault("sep", "\t")
    kwargs.setdefault("index", False)

    logger.info(f"Writing to {path}")

    df.to_csv(path, *args, **kwargs)

    return df


# The root logger acts as a sponge for loggers from external modules:
root_logger = setup_logger(name="")
# add_stream_handler(root_logger, logging.DEBUG)

# This is the main logger for the package
# Note that by default, it has no handlers. `setup_logger()` should be called after the
# `if __name__ == "__main__":` clause in each module.
logger = setup_logger()

log_step = pandas_utils.log_step(
    time_taken=False, shape_delta=False, print_fn=logger.info, display_args=False
)
