"""Module docstring."""

import logging
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

logger = logging.getLogger(__name__)


def main():
    """Run as script."""

    # Set plot style
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)
    palette = sns.color_palette()

    _ = pd.read_csv(C.GENE_LIST_GNOMAD_CST, sep="\t", header=None)

    pass


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
