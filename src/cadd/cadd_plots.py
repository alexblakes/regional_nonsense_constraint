"""Plot CADD score summary data."""

# Imports
import logging
import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from upsetplot import plot

import src
from src import constants as C
from src import visualisation as vis
from src.visualisation import maps_plots

# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_PALETTE = vis.color_palette("regions")[::-1]
_METRIC = "CADD Phred"
_PATHS = [C.STATS_CADD_SYN, C.STATS_CADD_MIS, C.STATS_CADD_NON]
_CSQS = ["Synonymous", "Missense", "Nonsense"]
_CONSTRAINT = ["constrained", "unconstrained"]
_FIGSIZE = (12 * C.CM, 6 * C.CM)

# Logging
logger = logging.getLogger(__name__)


# Functions
def read_data(path):
    return pd.read_csv(path, sep="\t", index_col="region")


def plot_cadd(x, y, ax=None, xerr=None, yerr=None, color=_PALETTE, **kwargs):
    # Default kwargs to ax.errorbar
    kwargs.setdefault("linestyle", "")

    # Get current axis if not specified
    if ax == None:
        ax = plt.gca()

    # Scatter plot
    ax.scatter(
        x=x,
        y=y,
        color=color,
    )

    # Error bars
    ax.errorbar(
        x=x,
        y=y,
        xerr=xerr,
        yerr=yerr,
        ecolor=color,
        **kwargs,
    )

    return None


def main():
    # Read data
    syn, mis, non = [read_data(p).assign(csq=c) for p, c in zip(_PATHS, _CSQS)]

    df = pd.concat([syn, mis, non])

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
    _ = list(itertools.product(_CONSTRAINT, _CSQS))
    logger.debug(_)

    for ax, (constraint, csq) in zip(axs, itertools.product(_CONSTRAINT, _CSQS)):
        m1 = df["constraint"] == constraint
        m2 = df["csq"] == csq
        data = df[m1 & m2]

        plot_cadd(data["mean"], data.index, xerr=data.ci_95, ax=ax)

    # s1 = syn[syn["constraint"] == "constrained"]

    # plot_cadd(x=s1["mean"], y=s1.index, ax=axs[0], xerr=s1["ci_95"])

    # axs[0].scatter(y=s1.index, x=s1["mean"])
    # axs[0].errorbar(
    #     y=s1.index,
    #     x=s1["mean"],
    #     xerr=s1["ci_95"],
    #     ecolor=_PALETTE,
    #     linestyle="",
    # )
    # axs[0].scatter(y=s1.index, x=s1["mean"], color=_PALETTE)
    # for ax, df in zip(axs[:len(_CSQS)], [syn, mis, non]):
    #     df = df[df["constraint"] == "constrained"]
    #     phylop_plots.horizontal_bars(df["mean"], ax=ax, xerr=df["ci_95"])

    # for ax, df in zip(axs[len(_CSQS):], [syn, mis, non]):
    #     df = df[df["constraint"] == "unconstrained"]
    #     phylop_plots.horizontal_bars(df["mean"], ax=ax, xerr=df["ci_95"])

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
    plt.savefig("data/plots/cadd.png", dpi=600)

    return syn, mis, non


if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)
    main()
