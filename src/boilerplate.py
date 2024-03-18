"""Boilerplate code for most modules."""

# Imports
import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C


# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"


# Logging
logger = logging.getLogger(__name__)


# Functions
def main():
    pass


if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)
    main()
