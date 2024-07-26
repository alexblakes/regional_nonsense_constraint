"""Create an upset plot of intersecting constrained regions."""

import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import seaborn as sns
import upsetplot

from src import constants as C

_FILE_IN = "data/statistics/constraint_upset_data.tsv"
_PNG = "data/plots/constraint/upset_plot.png"
_SVG = "data/plots/constraint/upset_plot.svg"


def read_upset_data(path):
    return pd.read_csv(path, sep="\t", index_col="enst")


def reformat_upset_data(df, indicators=C.REGION_LABELS):
    return upsetplot.from_indicators(indicators, data=df)


def upset_plot(upset_data, **kwargs):
    
    kwargs.setdefault("sort_by", "degree")
    kwargs.setdefault("sort_categories_by", "-input")
    kwargs.setdefault("max_subset_size", 2500)
    kwargs.setdefault("element_size", 18)
    kwargs.setdefault("intersection_plot_elements", 3)
    kwargs.setdefault("totals_plot_elements", 3)
    kwargs.setdefault("facecolor", "black")

    return upsetplot.UpSet(upset_data, **kwargs)


def customise_totals(ax=None):
    if not ax:
        ax = plt.gca()

    ax.grid(False)
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    ax.set_title("Total constrained", loc="right", ha="right")
    for c in ax.containers:
        ax.bar_label(c, padding=5)

    return ax


def customise_intersection(ax=None):
    if not ax:
        ax = plt.gca()

    ax.set_ylabel("Intersecting\nconstrained", va="bottom", loc="bottom")
    ax.grid(False)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    for c in ax.containers:
        ax.bar_label(c, padding=1, fontsize=7)

    return ax


def customise_color(upset):
    for region, colour in zip(C.REGION_LABELS, sns.color_palette()):
        upset.style_categories(
            region,
            shading_facecolor=colors.to_rgba(colour, alpha=0.5),
            bar_facecolor=colour,
        )
    return upset


def main():
    """Run as script."""

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])
    fig, ax = plt.subplots(1, 1, figsize=(18 * C.CM, 4.33 * C.CM), layout="constrained")
    ax.set_axis_off()

    upset_object = read_upset_data(_FILE_IN).pipe(reformat_upset_data).pipe(upset_plot)
    customise_color(upset_object)

    upset_axs = upset_object.plot(fig)

    customise_totals(upset_axs["totals"])
    customise_intersection(upset_axs["intersections"])

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return upset_axs


if __name__ == "__main__":
    main()
