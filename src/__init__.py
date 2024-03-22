""""Package docstring."""

import logging


# Formatters
_FORMATTER_ROOT = logging.Formatter(
    fmt="[%(asctime)s] %(levelname)s %(pathname)s %(funcName)s(): %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)
_FORMATTER_SRC = logging.Formatter(
    fmt="[%(asctime)s] %(levelname)s %(filename)s %(funcName)s(): %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)

# Log file paths
_LOGFILE_ROOT = "data/logs/_root.log"
_LOGFILE_SRC = "data/logs/_src.log"


def setup_logger(
    name="src",
    formatter=_FORMATTER_SRC,
    logfile=None,
    mode="w",
    stream=False,
    level=logging.DEBUG,
):
    """Create a new logger."""
    logger = logging.getLogger(name)

    # Log to file
    if logfile:
        file_handler = logging.FileHandler(logfile, mode=mode)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(level)
        logger.addHandler(file_handler)

    # Log to console
    if stream:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        stream_handler.setLevel(level)
        logger.addHandler(stream_handler)

    return logger


# Set up the root logger
root_logger = setup_logger(
    name="", logfile=_LOGFILE_ROOT, formatter=_FORMATTER_ROOT, stream=False, mode="a"
)
root_logger = setup_logger(
    name="", formatter=_FORMATTER_ROOT, stream=True, level=logging.INFO
)

# Set up the src logger
src_logger = setup_logger(logfile=_LOGFILE_SRC, mode="a")
