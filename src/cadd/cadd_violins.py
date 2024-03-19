"""Plot CADD score summary data."""

# Imports
import itertools
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis
from src.visualisation import maps_plots

# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_PALETTE = vis.color_palette("regions")[::-1]
_METRIC = "CADD Phred"
_CSQS = ["Synonymous", "Missense", "Nonsense"]
_CONSTRAINT = ["Constrained", "Unconstrained"]
_FIGSIZE = (12 * C.CM, 6 * C.CM)
_REGION_LABELS = ["Whole CDS", "NMD target", "Start proximal", "Long exon", "Distal"][
    ::-1
]  # Reversed for plotting
# Logging
logger = logging.getLogger(__name__)


# Functions
def read_data(path):
    return pd.read_csv(
        path,
        sep="\t",
        usecols=["csq", "cadd_phred", "region", "constraint"],
        dtype={"cadd_phred": np.float16},
        low_memory=False,
        # nrows=10000,
    )


def tidy_data(df):
    df = df.copy()

    # Constraint
    df = df.dropna(subset="constraint")
    df["constraint"] = df["constraint"].str.capitalize()

    # Consequences
    df["csq"] = df["csq"].replace(
        {
            "synonymous_variant": "Synonymous",
            "missense_variant": "Missense",
            "stop_gained": "Nonsense",
        }
    )

    # Region
    df["region"] = df["region"].replace(
        {
            "nmd_target": "NMD target",
            "start_proximal": "Start proximal",
            "long_exon": "Long exon",
            "distal_nmd": "Distal",
        }
    )

    return df


def get_whole_cds_data(df):
    cds = df.copy().assign(region="Whole CDS")
    return pd.concat([df, cds])


def figure(df):
    # Instantiate the figure
    fig, axs = plt.subplots(
        len(_CONSTRAINT),
        len(_CSQS),
        figsize=_FIGSIZE,
        layout="constrained",
        gridspec_kw={"hspace": 0.1, "wspace": 0.05},
        sharey=True,
        sharex="col",
    )
    axs = axs.flatten()

    # Violin plots
    for ax, (constraint, csq) in zip(axs, itertools.product(_CONSTRAINT, _CSQS)):
        
        # Split data
        m1 = df["constraint"] == constraint
        m2 = df["csq"] == csq
        data = df[m1 & m2]

        data_split = [
            data.loc[data["region"] == r, "cadd_phred"] for r in _REGION_LABELS
        ]

        # Plot violins
        violins = ax.violinplot(
            data_split,
            vert=False,
            widths=1,
            showmedians=False,
            showextrema=False,
        )

        # Customise violin color
        for body, color in zip(violins["bodies"], _PALETTE):
            body.set_facecolor(color)
            body.set_alpha(1)

        # Add box and whisker
        quantiles = [
            np.percentile(data_split[n], [25, 50, 75]) for n in range(len(data_split))
        ]
        q25 = [quantiles[n][0] for n in range(len(quantiles))]
        q50 = [quantiles[n][1] for n in range(len(quantiles))]
        q75 = [quantiles[n][2] for n in range(len(quantiles))]

        indices = np.arange(1, len(quantiles) + 1)

        ax.scatter(x=q50, y=indices, marker="o", color="white", s=2, zorder=3)
        ax.hlines(y=indices, xmin=q25, xmax=q75, color="black", linewidth=2)

        # Add y tick labels
        ax.set_yticks(
            [n + 1 for n in range(len(_REGION_LABELS))], labels=_REGION_LABELS
        )

    # Figure-level adjustments
        
    # Add x label
    for ax, csq in zip(axs[-len(_CSQS):], _CSQS):
        ax.set_xlabel(f"CADD Phred\n({csq})")

    # Add Axes titles
    for ax in axs[: len(_CSQS)]:
        ax.set_title("Constrained")
    
    for ax in axs[-len(_CSQS):]:
        ax.set_title("Unconstrained")

    # Trim x limits
    for ax in [axs[0], axs[3]]:
        ax.set_xlim(0, 20)

    for ax in [axs[1], axs[4]]:
        ax.set_xlim(10, 35)

    for ax in [axs[2], axs[5]]:
        ax.set_xlim(25, 50)


    # Save figure
    plt.savefig("data/plots/cadd.svg")
    plt.savefig("data/plots/cadd.png", dpi=600)
    
    plt.close()

    return None


def main():
    df = read_data(C.CADD_ANNOTATED).pipe(tidy_data).pipe(get_whole_cds_data)

    return df


if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)
    main()
