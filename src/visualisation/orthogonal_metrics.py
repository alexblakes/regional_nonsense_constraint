"""Plot phyloP, AlphaMissense, and pext scores in constrained and unconstrained 
regions.
"""

from collections import defaultdict
import logging
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/scratch/orthogonal_metrics_shuffled.tsv"
_PNG = "data/plots/orthogonal_metrics/phylop_alphamis_pext.png"
_SVG = "data/plots/orthogonal_metrics/phylop_alphamis_pext.svg"
_DTYPES = defaultdict(lambda: "float16", region="category", constraint="category")

logger = logging.getLogger(__name__)


def read_data(path: str) -> pd.DataFrame:
    logger.info("Reading data.")

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=_DTYPES,
        nrows=1000000,
    )

    return df


def sort_data(df: pd.DataFrame) -> pd.DataFrame:
    return df.assign(
        region=lambda x: pd.Categorical(
            x.region, categories=C.REGION_LABELS, ordered=True
        )
    ).sort_values("region")


def grouped_boxplot(df: pd.DataFrame, ax: plt.Axes = None, **kwargs):  # type: ignore
    if not ax:
        ax = plt.gca()

    kwargs.setdefault("x", "region")
    kwargs.setdefault("y", "phylop")
    kwargs.setdefault("hue", "constraint")
    kwargs.setdefault("showfliers", False)
    kwargs.setdefault("showmeans", True)
    kwargs.setdefault("whis", (5, 95))
    kwargs.setdefault("gap", 0.2)

    meanprops = dict(
        marker="o", markerfacecolor="black", markeredgecolor="black", markersize=3
    )

    sns.boxplot(
        df,
        ax=ax,
        meanprops=meanprops,
        palette=[
            sns.color_palette()[0],
            vis.adjust_lightness(sns.color_palette()[0], 1.25),
        ],
        linecolor="black",
        **kwargs,
    )

    ax.set_xlabel("")

    return ax


def customise_plot(ax=None, title=None, legend=False, **kwargs):
    if not ax:
        ax = plt.gca()

    ticks = ax.get_xticks()
    labels = ax.get_xticklabels()
    ax.set_xticks(ticks, labels, rotation=45, rotation_mode="anchor", ha="right")
    ax.tick_params(axis="x", length=0)

    ax.set(**kwargs)

    if legend:
        ax.legend(labelcolor="black", loc="lower left", bbox_to_anchor=(0, 1), ncols=1)
    else:
        ax.legend().remove()

    return ax


def recolor(ax: plt.Axes, palette=sns.color_palette()) -> plt.Axes:
    boxes = ax.findobj(patches.PathPatch)
    colors = list(palette) + [vis.adjust_alpha(p, 0.35) for p in palette]

    for box, color in zip(boxes, colors):
        box.set_facecolor(color)

    return ax


def main():
    """Run as script."""

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])
    fig, axs = plt.subplots(1, 3, figsize=(18 * C.CM, 6 * C.CM), layout="constrained")

    annotations = read_data(_FILE_IN).pipe(sort_data)

    metrics = ["phylop", "alpha_mis", "pext"]
    labels = ["phyloP", "AlphaMissense", "pext"]
    legends = [True, False, False]

    for ax, metric, label, legend in zip(axs, metrics, labels, legends):
        grouped_boxplot(annotations, ax, y=metric)
        customise_plot(ax, ylabel=label, legend=legend)
        recolor(ax)

    logger.info("Saving plot.")
    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return axs


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
