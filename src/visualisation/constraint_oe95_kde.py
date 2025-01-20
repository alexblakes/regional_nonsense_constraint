"""Create KDE plots for Nonsense OE95 per region."""



import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import statistics_for_plots as sfp

FILE_IN = "data/final/regional_nonsense_constraint.tsv"
PNG = "data/plots/constraint/oe95_kde.png"
SVG = "data/plots/constraint/oe95_kde.svg"

logger = src.logger


def parse_data(path=FILE_IN):
    return (
        pd.read_csv(FILE_IN, sep="\t")
        .loc[:, ["region", "oe_ci_hi"]]
        .dropna()
        .pipe(sfp.sort_column)
    )


def customise_plot(ax=None, title=None):
    ax = ax or plt.gca()

    ax.set_title(title, ma="center")
    # ax.text(1, 1, s=title, ma="left", ha="right", va="top", transform=ax.transAxes)
    ax.set_xlabel("Nonsense OE95")
    ax.set_ylabel(None)
    ax.set_yticks([])
    ax.set_xlim(0, 4)
    ax.axvline(0.6, ls="--", color="black")

    return ax


def annotate_plot(ax=None):
    ax = ax or plt.gca()

    # Annotate vertical line
    ax.text(s="0.6", x=0.55, y=0.95, ha="right", va="top", rotation=90)
    
    # Annotate constraint arrow
    ax.annotate(
        "",
        (0.8, 0.1),
        (4, 0.1),
        ax.transData,
        ax.transData,
        arrowprops=dict(arrowstyle="->,head_width=0.4"),
        ha="right",
        va="bottom",
        ma="right"
    )
    ax.text(
        s="Stronger\nconstraint",
        x=3.9,
        y=0.12,
        ha="right",
        va="bottom",
        ma="right",
        transform=ax.transData,
    )


def plot(axs):
    df = parse_data()

    groups = [g for _, g in df.groupby("region")]
    titles = df.region.unique()
    colors = sns.color_palette()

    for data, ax, title, color in zip(groups, axs, titles, colors):
        x = data["oe_ci_hi"]
        n = len(data)
        title = f"{title}\nN={n:,}"

        sns.kdeplot(ax=ax, x=x, gridsize=500, color=color, fill=True)
        customise_plot(ax, title)
    
    annotate_plot(axs[0])

    return axs


def main():
    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])
    fig, axs = plt.subplots(
        1, 5, figsize=(18 * C.CM, 4 * C.CM), layout="constrained", sharey=True
    )

    plot(axs)

    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    return None


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
