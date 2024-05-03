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
_FILE_IN_NMD_TARGET = "data/statistics/ora_nmd_target.tsv"
_FILE_IN_START_PROXIMAL = "data/statistics/ora_start_proximal.tsv"
_FILE_IN_LONG_EXON = "data/statistics/ora_long_exon.tsv"
_FILE_IN_DISTAL = "data/statistics/ora_distal.tsv"
_DROP_COLS = [
    "geneSet",
    "link",
    "size",
    "overlap",
    "expect",
    "pValue",
    "FDR",
    "overlapId",
    "userId",
]
_NMD_REGIONS_DICT = {
    "nmd_target": "NMD target",
    "start_proximal": "Start proximal",
    "long_exon":"Long exon",
    "distal":"Distal",
}
_SOURCES = ["hpo", "bp", "mf"]

logger = logging.getLogger(__name__)


def read_ora_data(paths):
    dfs = [pd.read_csv(path, sep="\t") for path in paths]
    return (
        pd.concat(dfs, axis=0)
        .drop(_DROP_COLS, axis=1)
        .replace({"region": _NMD_REGIONS_DICT})
    )


def rank_enrichment(df):
    g = df.groupby(["ref", "region", "db"])
    rank = lambda x: x.rank(method="first", ascending=False)
    df["enrichment_rank"] = g["enrichmentRatio"].transform(rank)

    return df.sort_values("enrichment_rank")


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
        width=data["enrichmentRatio"],
        tick_label=data["description"],
        color=color,
    )

    ax.bar_label(b, fmt="{:.2f}", padding=3)
    ax.invert_yaxis()
    ax.spines[["left"]].set_visible(False)
    ax.set_xlabel("Fold-enrichment")
    ax.label_outer()

    # Shorten very long tick labels
    shorten = lambda x: textwrap.shorten(x.get_text(), 60)
    ax.set_yticklabels(map(shorten, ax.get_yticklabels()))

    return ax


def plot_enrichment(data, bg_name, queries, nrows=4, palette=sns.color_palette()):
    """Create gene set enrichment figures."""

    grouped_data = data.groupby(["db", "region"])

    for source in _SOURCES:
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

        plt.savefig(f"data/plots/gene_set_enrichment/{bg_name}_{source}.png", dpi=600)
        plt.savefig(f"data/plots/gene_set_enrichment/{bg_name}_{source}.svg")


def main():
    """Run as script."""

    # Plotting styles
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)

    # Parse ORA data
    paths = (
        _FILE_IN_NMD_TARGET,
        _FILE_IN_START_PROXIMAL,
        _FILE_IN_LONG_EXON,
        _FILE_IN_DISTAL,
    )
    df = read_ora_data(paths).pipe(rank_enrichment)

    # Stratify by "background" gene set
    group_bg = df.groupby("ref")
    df_all = group_bg.get_group("all").copy()
    df_gnomad = group_bg.get_group("gnomad").copy()

    # Plot enrichment
    dfs = [df_all, df_gnomad]
    names = ["all", "gnomad"]
    query = C.NMD_REGION_LABELS
    nrows = len(query)
    palette = sns.color_palette()[1:]

    for (
        df,
        name,
    ) in zip(dfs, names):
        plot_enrichment(df, name, query, nrows, palette)

    plt.close("all")

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
