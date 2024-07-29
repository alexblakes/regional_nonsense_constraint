"""Plot phyloP, AlphaMissense, and pext scores in constrained and unconstrained 
regions.
"""

from collections import defaultdict
import logging
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import patches
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/statistics/orthogonal_metrics.tsv.gz"
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
    kwargs.setdefault("palette", ["0.7", "0.9"])
    kwargs.setdefault("linecolor", "black")

    meanprops = dict(
        marker="o", markerfacecolor="black", markeredgecolor="black", markersize=3
    )

    sns.boxplot(df, ax=ax, meanprops=meanprops, **kwargs)

    ax.set_xlabel("")
    ax.legend().remove()

    return ax


def customise_plot(ax: plt.Axes = None, legend: bool = False, **kwargs):
    if not ax:
        ax = plt.gca()

    ax.set(**kwargs)

    ax.set_xticks(
        ticks=ax.get_xticks(),
        labels=ax.get_xticklabels(),
        rotation=45,
        rotation_mode="anchor",
        ha="right",
    )

    if legend:
        ax.legend(labelcolor="black", loc="lower left", bbox_to_anchor=(0, 1), ncols=1)

    return ax


def recolor(ax: plt.Axes, palette) -> plt.Axes:
    boxes = ax.findobj(patches.PathPatch)

    dark_colors = list(palette)
    light_colors = [vis.adjust_alpha(p, 0.35) for p in palette]
    colors = dark_colors + light_colors

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
        recolor(ax, sns.color_palette())

    logger.info("Saving plot.")
    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return axs


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
