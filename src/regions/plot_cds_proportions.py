"""Plot ACMG classification of nonsense / frameshift variants in ClinVar."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import pandas as pd
import seaborn as sns

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/statistics/regions_cds_proportions.tsv"
_PNG = "data/plots/regions_cds_proportions.png"
_SVG = "data/plots/regions_cds_proportions.svg"

logger = logging.getLogger(__name__)


def read_cds_proportions(path):
    return pd.read_csv(path, sep="\t", index_col="region")


def plot_bars(series, ax=None, **kwargs):
    kwargs.setdefault("color", sns.color_palette())

    if not ax:
        ax = plt.gca()

    x = series.index

    ax.bar(x=x, height=series, **kwargs)

    return ax


def customise_plot(ax=None):
    if not ax:
        ax = plt.gca()

    # x axis
    labels = ax.get_xticklabels()
    ticks = np.arange(len(labels))
    ax.set_xticks(ticks, labels, rotation=45, ha="right", rotation_mode="anchor")

    # y axis
    ax.set_ylabel("Proportion of CDS")
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))

    # Bar labels
    bars = ax.containers[0]
    ax.bar_label(bars, fmt="{:.0%}")

    return ax


def main():
    """Run as script."""
    cds_proportions = read_cds_proportions(_FILE_IN)

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])

    fig, ax = plt.subplots(1, 1, figsize=(5 * C.CM, 5 * C.CM), layout="constrained")

    plot_bars(cds_proportions["proportion"])
    customise_plot()

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)

    return ax


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
