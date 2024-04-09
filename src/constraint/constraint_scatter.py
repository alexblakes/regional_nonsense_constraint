"""Module docstring."""

import itertools
import logging
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
        layout="constrained",
        # sharex="col",
        # sharey="row",
        subplot_kw=dict(box_aspect=1),
    )
    axs = axs.flatten("F")

    # Create plots
    for ax, (xlabel, ylabel) in zip(axs, subsets):
        if xlabel == ylabel:
            continue

        # Keep only relevant columns
        x_constraint = xlabel + "_constraint"
        y_constraint = ylabel + "_constraint"
        x_oe = xlabel + "_oe"
        y_oe = ylabel + "_oe"

        data = df[[x_constraint, x_oe, y_constraint, y_oe]].copy()

        # Keep entries where at least one region is constrained
        m1 = data[x_constraint] == "constrained"
        m2 = data[y_constraint] == "constrained"

        data = data[m1 | m2]

        # Drop entries with NaNs in OE
        data = data.dropna(subset=[y_oe])

        # Plot the data
        x = data[x_oe]
        y = data[y_oe]
        ax.scatter(x, y)

        ax.set_xlabel(region_dict[xlabel])
        ax.set_ylabel(region_dict[ylabel])
        # ax.axline(xy1=(0,0), slope=1)
        # ax.set_aspect(np.diff(ax.get_xlim()) / np.diff(ax.get_ylim()))
        # ax.label_outer()

        # break

    plt.savefig("data/plots/_", dpi=600)
    # plt.savefig("data/plots/_.svg")
    plt.close("all")

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
