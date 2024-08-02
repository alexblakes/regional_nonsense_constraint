"""Boilerplate code for most modules."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).relative_to(Path.cwd()).with_suffix('.log').parts)}"

logger = logging.getLogger(__name__)


def main():
    """Run as script."""
    return None


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
