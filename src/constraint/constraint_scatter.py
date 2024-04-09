"""Module docstring."""

import itertools
import logging
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import adjustText

import src
from src import constants as C
from src.constraint import gene_lists

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

logger = logging.getLogger(__name__)


def main():
    """Run as script."""

    # Read data
    df = pd.read_csv(
        C.REGIONAL_NONSENSE_CONSTRAINT,
        sep="\t",
        usecols=[
            "enst",
            "region",
            "n_exp",
            "oe",
            "z",
            "pli",
            "loeuf",
            "syn_z",
            "constraint",
        ],
    )

    # Filter data
    df = df.dropna(subset="z")
    df = df[df["syn_z"] > -1]
    df = df[df["n_exp"] >= 5]

    # Annotate with gene symbols
    gene_ids = pd.read_csv(
        C.CANONICAL_CDS_GENE_IDS,
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
            values=["constraint", "oe"],
        )
        .reorder_levels([1, 0], axis=1)
        .sort_index(axis=1)
    )
    df.columns = ["_".join([x, y]) for x, y in df.columns]

    # Set plot style
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)

    # Get parameters for plots
    _REGION_NAMES = ["nmd_target", "start_proximal", "long_exon", "distal_nmd"]
    _REGION_LABELS = ["NMD target", "Start proximal", "Long exon", "Distal"]
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
        if xlabel == ylabel:
            continue

        # Keep only relevant columns
        x_constraint = xlabel + "_constraint"
        y_constraint = ylabel + "_constraint"
        x_oe = xlabel + "_oe"
        y_oe = ylabel + "_oe"

        data = df[[x_constraint, x_oe, y_constraint, y_oe]].copy().reset_index("symbol")

        # Keep entries where at least one region is constrained
        m1 = data[x_constraint] == "constrained"
        m2 = data[y_constraint] == "constrained"

        data = data[m1 | m2]

        # Drop entries with NaNs in OE
        data = data.dropna(subset=[y_oe])

        # Plot scatter plots
        x = data[x_oe]
        y = data[y_oe]
        ax.scatter(x, y, alpha=0.5, linewidth=0)

        ax.set_xlabel(region_dict[xlabel])
        ax.set_ylabel(region_dict[ylabel])
        # ax.axline(xy1=(0,0), slope=1)
        # ax.set_aspect(np.diff(ax.get_xlim()) / np.diff(ax.get_ylim()))
        ax.label_outer()

        # Highlight the strongest outliers
        data["diff"] = (data[x_oe] - data[y_oe]).astype(float)
        # data = data.dropna(subset="diff")
        min3 = data.nsmallest(3, "diff")
        max3 = data.nlargest(3, "diff")

        print(max3)

        ax.scatter(
            min3[x_oe], min3[y_oe], color="None", edgecolor="black", linewidth=0.5
        )
        ax.scatter(
            max3[x_oe], max3[y_oe], color="None", edgecolor="black", linewidth=0.5
        )

        annots = []
        for i, row in min3.iterrows():
            annots.append(
                ax.text(
                    x=row[x_oe],
                    y=row[y_oe],
                    s=row["symbol"],
                    ha="left",
                    va="center",
                    size=7
                )
            )
        adjustText.adjust_text(
            annots,
            ax=ax,
            expand=(1.2, 1.5),
            avoid_self=True,
            arrowprops=dict(arrowstyle="-", linewidth=0.5),
        )

        # break

    plt.savefig("data/plots/_", dpi=600)
    # plt.savefig("data/plots/_.svg")
    plt.close("all")

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
