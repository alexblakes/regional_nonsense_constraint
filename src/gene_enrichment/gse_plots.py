"""Plot gene set enrichment comparisons."""

import itertools
import logging
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_SOURCES = ["HP", "GO:BP", "GO:MF"]
_SOURCE_NAMES = ["hpo", "gobp", "gomf"]
_ALL_QUERIES = ["gnomAD", "NMD target", "Start proximal", "Long exon", "Distal"]
_GNOMAD_QUERIES = ["NMD target", "Start proximal", "Long exon", "Distal"]

logger = logging.getLogger(__name__)


def create_figure(
    nrows=5,
    ncols=1,
    figsize=(12 * C.CM, 18 * C.CM),
    sharex=True,
    layout="constrained",
    **kwargs,
):
    """Instantiate a figure with multiple subplots."""

    fig, axs = plt.subplots(
        nrows,
        ncols,
        figsize=figsize,
        sharex=sharex,
        layout=layout,
        **kwargs,
    )
    axs = axs.flatten()

    return fig, axs


def enrichment_bars(ax, data, color):
    """Plot enrichment bars for a given data slice."""

    b = ax.barh(
        y=data["enrichment_rank"],
        width=data["enrichment"],
        tick_label=data["name"],
        color=color,
    )

    ax.bar_label(b, fmt="{:.2f}", padding=3)
    ax.invert_yaxis()
    ax.spines[["left"]].set_visible(False)
    ax.set_xlabel("Fold-enrichment")

    # Shorten very long tick labels
    shorten = lambda x: textwrap.shorten(x.get_text(), 60)
    ax.set_yticklabels(map(shorten, ax.get_yticklabels()))

    return ax


def plot_enrichment(data, bg_name, queries, nrows=5, palette=sns.color_palette()):
    """Create gene set enrichment figures."""

    grouped_data = data.groupby(["source", "query"])

    for source, source_name in zip(_SOURCES, _SOURCE_NAMES):
        fig, axs = create_figure(nrows=nrows)
        palette = itertools.cycle(palette)

        for query, ax in zip(queries, axs):
            color = next(palette)
            ax.set_title(query, loc="left", color=color)
            ax.set_yticks([])

            try:
                data_subset = grouped_data.get_group((source, query)).copy()
                enrichment_bars(ax, data_subset, color)

            except:
                continue

        plt.savefig(
            f"data/plots/gene_set_enrichment/{bg_name}_{source_name}.png", dpi=600
        )
        plt.savefig(f"data/plots/gene_set_enrichment/{bg_name}_{source_name}.svg")


def main():
    """Run as script."""

    # Plotting styles
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)

    # Read data
    df = pd.read_csv(C.STATS_GENE_SET_ENRICHMENT, sep="\t")

    # Stratify by "background" gene set
    group_bg = df.groupby("background")
    df_all = group_bg.get_group("All genes").copy()
    df_gnomad = group_bg.get_group("gnomAD").copy()

    # Plot enrichment
    dfs = [df_all, df_gnomad]
    names = ["all", "gnomad"]
    queries = [_ALL_QUERIES, _GNOMAD_QUERIES]
    nrows = [5, 4]
    palettes = [sns.color_palette(), sns.color_palette()[1:]]

    for df, name, query, _nrows, palette in zip(dfs, names, queries, nrows, palettes):
        plot_enrichment(df, name, query, _nrows, palette)

    plt.close("all")

    return None


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
