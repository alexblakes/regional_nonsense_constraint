""""Package docstring."""

import logging

from src import constants as C

# Module constants
_FORMATTER = logging.Formatter(
    fmt="[%(asctime)s] %(levelname)s %(name)s %(funcName)s(): %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)


# Functions
def add_filehandler(logger, logfile, level=logging.DEBUG):
    
    file_handler = logging.FileHandler(logfile)
    file_handler.setLevel(level)
    file_handler.setFormatter(_FORMATTER)
    logger.addHandler(file_handler)

    return logger

def add_streamhandler(logger, level=logging.INFO):
    
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(_FORMATTER)
    logger.addHandler(stream_handler)
    
    return logger

def setup_logger(
    name,
    propagate=True,
    logfile=None,
):
    """Create a new logger."""

    logger = logging.getLogger(name)
    logger.propagate = propagate
    logger.setLevel(logging.DEBUG)

    # Log to file
    if logfile:
        logger = add_filehandler(logger, logfile)

    # Log to std err
    logger = add_streamhandler(logger)

    return logger


root_logger = setup_logger("", logfile="data/logs/_all_python_modules.log")
