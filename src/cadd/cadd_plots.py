"""Plot CADD score summary data."""

# Imports
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import src
from src import constants as C
from src import visualisation as vis
from src.visualisation import phylop_plots

# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_PALETTE = vis.color_palette("maps")
_METRIC = "CADD Phred"
_CSQS = ["Synonymous", "Missense", "Nonsense"]
_CONSTRAINT = ["constrained", "unconstrained"]
_FIGSIZE = (12 * C.CM, 6 * C.CM)

# Logging
logger = logging.getLogger(__name__)


# Functions
def read_data(path):
    return pd.read_csv(path, sep="\t", index_col="region")


def main():
    # Read data
    syn = read_data(C.STATS_CADD_SYN)
    mis = read_data(C.STATS_CADD_MIS)
    non = read_data(C.STATS_CADD_NON)

    # Instantiate the figure
    fig, axs = plt.subplots(
        len(_CONSTRAINT),
        len(_CSQS),
        figsize=_FIGSIZE,
        layout="constrained",
        gridspec_kw={"hspace": 0.1, "wspace": 0.05},
        sharey=True,
        sharex="col",
    )
    axs = axs.flatten()

    # Plots
    for ax, df in zip(axs[:len(_CSQS)], [syn, mis, non]):
        df = df[df["constraint"] == "constrained"]
        phylop_plots.horizontal_bars(df["mean"], ax=ax, xerr=df["ci_95"])
    
    for ax, df in zip(axs[len(_CSQS):], [syn, mis, non]):
        df = df[df["constraint"] == "unconstrained"]
        phylop_plots.horizontal_bars(df["mean"], ax=ax, xerr=df["ci_95"])

    # # X axis labels in lowest Axes
    # for ax in axs[-len(_CSQS) :]:
    #     ax.set_xlabel(_METRIC)

    # # Constraint annotations
    # for ax, title in zip(axs, [x.capitalize() for x in sorted(_CONSTRAINT * len(_METRICS))]):
    #     ax.set_title(title)

    # # Significance annotations
    # for ax in axs[:3]:
    #     bars = ax.containers[1]
    #     ax.bar_label(bars, labels=[r"$\star$"] * len(bars))
    #     # for bar in ax.containers[1]:
    #     #     ax.bar_label(bar, labels="*")

    # Save figure
    plt.savefig("data/plots/cadd.svg")
    plt.savefig("data/plots/cadd.png", dpi = 600)

    return syn, mis, non


if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)
    main()
