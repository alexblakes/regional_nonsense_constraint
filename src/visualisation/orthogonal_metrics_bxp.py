"""Plot phyloP, AlphaMissense, and pext scores in constrained and unconstrained 
regions.
"""



import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import patches
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis

_FILE_IN = "data/statistics/orthogonal_metrics_bxp_stats.tsv"
_PNG = "data/plots/orthogonal_metrics/boxplots.png"
_SVG = "data/plots/orthogonal_metrics/boxplots.svg"

logger = src.logger


def read_data(path: str = _FILE_IN) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def customise_plot(ax=None, legend=False, legend_kwargs={}, **kwargs):
    legend_kwargs.setdefault("loc", "lower left")
    legend_kwargs.setdefault("bbox_to_anchor", (0, 1))
    legend_kwargs.setdefault("labelcolor", "black")

    if not ax:
        ax = plt.gca()

    ax.set_xlim(0)

    ax.set_xticks(
        ticks=ax.get_xticks(),
        labels=ax.get_xticklabels(),
        rotation=45,
        rotation_mode="anchor",
        ha="right",
    )

    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6, min_n_ticks=3))

    # Add legend
    if legend:
        color = sns.color_palette()[0]
        light_color = vis.adjust_alpha(color, 0.35)
        patch1 = patches.Patch(
            facecolor=color, edgecolor="black", linewidth=0.5, label="Constrained"
        )
        patch2 = patches.Patch(
            facecolor=light_color,
            edgecolor="black",
            linewidth=0.5,
            label="Unconstrained",
        )

        ax.legend(handles=[patch2, patch1], **legend_kwargs)

    ax.set(**kwargs)

    return ax


def recolor(ax: plt.Axes, palette) -> plt.Axes:
    boxes = ax.findobj(patches.PathPatch)

    dark = list(palette)
    light = [vis.adjust_alpha(p, 0.35) for p in palette]
    colors = dark + light

    for box, color in zip(boxes, colors):
        box.set_facecolor(color)

    return ax

def plot(axs):
    
    df = read_data().set_index(["metric", "region", "constraint"])

    subsets = [g for _, g in df.groupby(level="metric", sort=False)]
    legends = [1, 0, 0, 0]

    for subset, ax, legend in zip(subsets, axs, legends):
        groups = [g for _, g in subset.groupby("constraint")]
        start_positions = [2, 1]

        for data, start in zip(groups, start_positions):
            data = data.to_dict(orient="records")
            n_boxes = len(data)
            positions = [start + i * 3 for i in range(n_boxes)]

            ax.bxp(
                data,
                positions,
                showfliers=False,
                showmeans=True,
                widths=0.8,
                meanprops=dict(marker="o", markerfacecolor="black", markersize=3),
                boxprops=dict(lw=0.5),
                medianprops=dict(color="black", lw=0.5),
                whiskerprops=dict(lw=0.5),
                capprops=dict(lw=0.5),
                patch_artist=True,
            )

        ylabel = subset.index.get_level_values("metric")[0]
        xticks = [1.5 + i * 3 for i in range(len(C.REGIONS))]
        xticklabels = subset.index.get_level_values("region").unique()

        customise_plot(
            ax,
            legend,
            xticks=xticks,
            xticklabels=xticklabels,
            ylabel=ylabel,
            legend_kwargs=dict(ncols=2)
        )
        recolor(ax, sns.color_palette())

        asterisk_xs = [i for i in range(16) if not i % 3 == 0]
        asterisk_ys = subset["whishi"]
        asterisk_ps = [0,1] * 5
        vis.add_significance_asterisk(asterisk_xs, asterisk_ys, asterisk_ps, ax=ax, y_adj=4)
    
    return None


def main():
    """Run as script."""

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])
    fig, axs = plt.subplots(2, 2, figsize=(12 * C.CM, 12 * C.CM), layout="constrained")
    axs = axs.flatten()

    plot(axs)

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return fig, axs


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
