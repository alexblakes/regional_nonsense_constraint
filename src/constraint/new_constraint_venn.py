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

    gene_set = lambda x: set(pd.read_csv(x, header=None).iloc[:,0].to_list())

    gnomad = gene_set(C.GENE_LIST_GNOMAD_CST)
    target = gene_set(C.GENE_LIST_NMD_TARGET)
    start = gene_set(C.GENE_LIST_START_PROX)
    long_exon = gene_set(C.GENE_LIST_LONG_EXON)
    distal = gene_set(C.GENE_LIST_DISTAL)

    return gnomad


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
