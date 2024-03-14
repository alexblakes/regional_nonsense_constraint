""""Package docstring."""

import logging

from src import constants as C

# Module constants
_FORMATTER_ROOT = logging.Formatter(
    fmt="[%(asctime)s] %(levelname)s %(pathname)s %(funcName)s(): %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)
_FORMATTER_SRC = logging.Formatter(
    fmt="[%(asctime)s] %(levelname)s %(filename)s %(funcName)s(): %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)
_LOGFILE_ROOT = "data/logs/_root.log"
_LOGFILE_SRC = "data/logs/_src.log"


# Functions
def add_filehandler(
    logger, logfile, level=logging.DEBUG, mode="a", formatter=_FORMATTER_SRC
):
    file_handler = logging.FileHandler(logfile, mode=mode)
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)

    logger.addHandler(file_handler)

    return None


def add_streamhandler(logger, level=logging.INFO, formatter=_FORMATTER_SRC):
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return None


def setup_logger(
    name,
    propagate=True,
    stream=False,
    logfile=None,
    file_formatter=_FORMATTER_SRC,
    stream_formatter=_FORMATTER_SRC,
):
    """Create a new logger."""

    logger = logging.getLogger(name)
    logger.propagate = propagate
    logger.setLevel(logging.DEBUG)

    # Log to file
    if logfile:
        add_filehandler(logger, logfile, formatter=file_formatter)

    # Log to std err
    if stream:
        add_streamhandler(logger, formatter=stream_formatter)

    return logger


def module_logger(logfile, name="src"):
    logger = logging.getLogger(name)
    add_filehandler(logger, logfile, mode="w")

    return logger

# Set up loggers
root_logger = setup_logger(
    "",
    logfile=_LOGFILE_ROOT,
    stream=True,
    file_formatter=_FORMATTER_ROOT,
    stream_formatter=_FORMATTER_SRC,
)
src_logger = setup_logger("src", logfile=_LOGFILE_SRC, stream=False)
