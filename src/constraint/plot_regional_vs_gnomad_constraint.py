"""Explore the relationship between nonsense O/E values and LOEUF scores."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/final/regional_nonsense_constraint.tsv"
OE_SVG = "data/plots/constraint/oe_vs_loeuf.svg"
OE_PNG = "data/plots/constraint/oe_vs_loeuf.png"
OE_CI_SVG = "data/plots/constraint/oe_ci_hi_vs_loeuf.svg"
OE_CI_PNG = "data/plots/constraint/oe_ci_hi_vs_loeuf.png"


logger = logging.getLogger(__name__)


def main():
    """Run as script."""

    df = pd.read_csv(_FILE_IN, sep="\t")

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
