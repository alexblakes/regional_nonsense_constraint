""""Package docstring."""

import logging

from src import constants as C


# Functions
def setup_logger(name, propagate=True, logfile=None,):
    """Create a new logger."""
    
    logger = logging.getLogger(name)
    logger.propagate = propagate
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter(
        fmt="[%(asctime)s] %(levelname)s %(name)s %(funcName)s(): %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
    )

    # Log to file
    if logfile:
        file_handler = logging.FileHandler(logfile)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # Log to std err
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger

root_logger = setup_logger("", logfile="data/logs/_all_python_modules.log")

