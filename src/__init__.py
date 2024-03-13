""""Package docstring."""

import logging

from src import constants as C

# Module constants
_FORMATTER = logging.Formatter(
    fmt="[%(asctime)s] %(levelname)s %(filename)s %(funcName)s(): %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)
_LOGFILE_ROOT = "data/logs/_root.log"
_LOGFILE_SRC = "data/logs/_src.log"


# Functions
def add_filehandler(logger, logfile, level=logging.DEBUG, mode="a"):
    file_handler = logging.FileHandler(logfile, mode=mode)
    file_handler.setLevel(level)
    file_handler.setFormatter(_FORMATTER)

    logger.addHandler(file_handler)

    return None


def add_streamhandler(logger, level=logging.INFO):
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(_FORMATTER)
    logger.addHandler(stream_handler)

    return None


def setup_logger(name, propagate=True, stream=False, logfile=None):
    """Create a new logger."""

    logger = logging.getLogger(name)
    logger.propagate = propagate
    logger.setLevel(logging.DEBUG)

    # Log to file
    if logfile:
        add_filehandler(logger, logfile)

    # Log to std err
    if stream:
        add_streamhandler(logger)

    return logger


def module_logger(logfile, name="src"):
    logger = logging.getLogger(name)
    add_filehandler(logger, logfile, mode="w")

    return logger


root_logger = setup_logger("", logfile=_LOGFILE_ROOT, stream=True)
src_logger = setup_logger("src", logfile=_LOGFILE_SRC, stream=False)
