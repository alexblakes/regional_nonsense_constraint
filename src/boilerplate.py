"""Boilerplate code for most modules."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"

logger = logging.getLogger(__name__)


def main():
    """Run as script."""
    pass


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
