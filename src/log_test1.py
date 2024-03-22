"""Test logging."""

# Imports
import logging
from pathlib import Path

import src
from src import log_test2

# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"


# Logging
logger = logging.getLogger(__name__)


# Functions
def main():
    """Run as script."""
    logger.info("Speaking from module 1.")
    log_test2.speak("log_test2 speaking from module 1.")


if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)
    main()
