"""Plot ACMG classification of nonsense / frameshift variants in ClinVar."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import pandas as pd

import src
from src import constants as C
from src.visualisation import bars

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/statistics/clinvar_vus_by_region.tsv"
_PNG = "data/plots/clinvar/clinvar_vus_by_region.png"
_SVG = "data/plots/clinvar/clinvar_vus_by_region.svg"

logger = logging.getLogger(__name__)


def read_vus_proportions(path):
    return pd.read_csv(path, sep="\t", index_col="region")


def customise_plot(ax=None):
    if not ax:
        ax = plt.gca()

    # x axis
    labels = ax.get_xticklabels()
    ticks = np.arange(len(labels))
    ax.set_xticks(ticks, labels, rotation=45, ha="right", rotation_mode="anchor")

    # y axis
    ax.set_ylabel("Proportion of PTVs\nwhich are VUS in ClinVar")
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))

    # Bar labels
    bars = ax.containers[0]
    ax.bar_label(bars, fmt="{:.0%}")

    return ax


def main():
    """Run as script."""
    vus_proportions = read_vus_proportions(_FILE_IN)

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])

    fig, ax = plt.subplots(1, 1, figsize=(5 * C.CM, 5 * C.CM), layout="constrained")

    bars.vertical_bars(vus_proportions["proportion_vus"])
    customise_plot()

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)

    return ax


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
