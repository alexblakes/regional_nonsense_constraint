"""Create a scatter plot for transcript-level O/E statistics."""

import logging

import matplotlib.pyplot as plt
from matplotlib import ticker
import pandas as pd

import src
from src import constants as C
from src import statistics_for_plots as sfp

FILE_IN = "data/final/regional_constraint_stats.tsv"
CSQS = {
    "synonymous_variant": "Synonymous",
    "missense_variant": "Missense",
    "stop_gained": "Nonsense",
}
PNG = "data/plots/constraint/oe_scatter.png"
SVG = "data/plots/constraint/oe_scatter.svg"

logger = logging.getLogger(__name__)


def parse_data(path=FILE_IN):
    return (
        pd.read_csv(FILE_IN, sep="\t")
        .query("region == 'transcript'")
        .query("enst != 'ENST00000589042'")  # drop TTN
        .loc[:, ["csq", "n_obs", "n_exp"]]
        .pipe(sfp.sort_column, column="csq", labels=CSQS)
    )


def customise_plot(ax=None, title=None):
    ax = ax or plt.gca()

    ax.axline((0, 0), slope=1, ls="--", color="grey", zorder=-1)
    ax.axis("square")
    ax.set_title(title)
    ax.set_xlabel("Expected")
    ax.set_ylabel("Observed")
    ax.xaxis.set_major_locator(ticker.MaxNLocator(3, min_n_ticks=3))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(3, min_n_ticks=3))

    return ax


def plot(fig, axs):
    # Get data
    df = parse_data()

    groups = [g for _, g in df.groupby("csq")]
    titles = df.csq.unique()

    # Create hexbin plots
    for ax, data, title in zip(axs, groups, titles):
        x = data["n_exp"]
        y = data["n_obs"]

        ax.hexbin(x, y, bins="log", cmap="Greys_r", linewidths=0, gridsize=40)
        customise_plot(ax, title)

    # Add colorbar
    hb = axs[-1].collections[-1]  # Get mappable for colorbar
    fig.colorbar(
        hb,
        ax=axs,
        label="Transcripts",
        shrink=0.6,
        ticks=ticker.LogLocator(numticks=10),
    )


def main():
    plt.style.use(C.STYLE_DEFAULT)
    fig, axs = plt.subplots(
        1,
        3,
        figsize=(12 * C.CM, 4 * C.CM),
        layout="constrained",
    )

    plot(fig, axs)

    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    return None


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
