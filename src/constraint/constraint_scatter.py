"""Plot constraint pairplots."""

import itertools
import logging
from pathlib import Path

import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import adjustText

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_REGIONAL_NONSENSE_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
_GENE_IDS = "data/interim/gene_ids.tsv"
_PNG = "data/plots/constraint/region_pair_plots.png"
_SVG = "data/plots/constraint/region_pair_plots.svg"
_REGION_NAMES = ["nmd_target", "start_proximal", "long_exon", "distal_nmd"]
_REGION_LABELS = ["NMD target", "Start proximal", "Long exon", "Distal"]

logger = logging.getLogger(__name__)


def main():
    """Run as script."""

    # Read data
    df = pd.read_csv(
        _REGIONAL_NONSENSE_CONSTRAINT,
        sep="\t",
        usecols=[
            "enst",
            "region",
            "n_exp",
            "oe",
            "oe_ci_hi",
            "p",
            "pli",
            "loeuf",
            "syn_p",
            "constraint",
        ],
    )

    # Filter data
    df = df.dropna(subset="p")
    df = df[df["syn_p"] >= stats.norm.cdf(-1)]
    df = df[df["n_exp"] >= 5]

    # Annotate with gene symbols
    gene_ids = pd.read_csv(
        _GENE_IDS,
        sep="\t",
        usecols=["transcript_id", "gene_name"],
        header=0,
    ).set_axis(["enst", "symbol"], axis=1)

    df = df.merge(gene_ids, validate="many_to_one")

    # Pivot table
    df = (
        df.pivot(
            index=["enst", "symbol"],
            columns="region",
            values=["constraint", "oe_ci_hi"],
        )
        .reorder_levels([1, 0], axis=1)
        .sort_index(axis=1)
    )
    df.columns = ["_".join([x, y]) for x, y in df.columns]

    # Set plot style
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)

    # Get parameters for plots
    region_dict = {a: b for a, b in zip(_REGION_NAMES, _REGION_LABELS)}
    dim = len(_REGION_NAMES)
    subsets = itertools.product(_REGION_NAMES, _REGION_NAMES)

    # Instantiate the figure
    fig, axs = plt.subplots(
        dim,
        dim,
        figsize=(18 * C.CM, 18 * C.CM),
        layout="tight",
        sharex="col",
        sharey="row",
        subplot_kw=dict(box_aspect=1),
    )
    axs = axs.ravel("F")

    # Create plots
    for ax, (xlabel, ylabel) in zip(axs, subsets):
        ax.set_xlabel(f"O/E\n{region_dict[xlabel]}")
        ax.set_ylabel(f"O/E\n{region_dict[ylabel]}")
        ax.label_outer()

        if xlabel == ylabel:
            continue

        # Keep only relevant columns
        x_constraint = xlabel + "_constraint"
        y_constraint = ylabel + "_constraint"
        x_oe = xlabel + "_oe_ci_hi"
        y_oe = ylabel + "_oe_ci_hi"

        data = df[[x_constraint, x_oe, y_constraint, y_oe]].copy().reset_index("symbol")

        # Keep entries where at least one region is constrained
        m1 = data[x_constraint] == "constrained"
        m2 = data[y_constraint] == "constrained"

        data = data[m1 | m2]

        # Plot scatter plots
        x = data[x_oe]
        y = data[y_oe]
        ax.scatter(x, y, alpha=0.5, linewidth=0)

        # Highlight the strongest outliers
        data["diff"] = (data[x_oe] - data[y_oe]).astype(float)
        min3 = data.nsmallest(5, "diff")

        ax.scatter(
            min3[x_oe], min3[y_oe], color="None", edgecolor="black", linewidth=0.5
        )

        # Annotate the N largest outliers
        annots = []
        for i, row in min3.iterrows():
            annots.append(
                ax.text(
                    x=row[x_oe],
                    y=row[y_oe],
                    s=row["symbol"],
                    ha="left",
                    va="top",
                    size=7,
                )
            )
        adjustText.adjust_text(
            annots,
            ax=ax,
            expand=(1.1, 1.2),
            time_lim=3,
            expand_axes=False,
            only_move="x+y",
            avoid_self=True,
            arrowprops=dict(arrowstyle="-", linewidth=0.5),
        )

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
