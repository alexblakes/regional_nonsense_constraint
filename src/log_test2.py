"""Test logging."""

# Imports
import logging
from pathlib import Path

import src

# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"


# Logging
logger = logging.getLogger(__name__)


# Functions
def speak(text):
    logger.info(text)

def main():
    """Run as script."""
    logger.info("Speaking from module 2.")


if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)
    main()
