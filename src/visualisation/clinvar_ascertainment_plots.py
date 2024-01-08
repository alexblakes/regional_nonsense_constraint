"""Docstring."""

# Imports
from pathlib import Path

import numpy as np
import pandas as pd

from src import constants as C
from src import setup_logger
from src.visualisation import color_tools as ct


# Logging
logger = setup_logger(Path(__file__).stem)


# Color palettes
VIBRANT_PALETTE = ct.color_palette("default")  # Vibrant palette
REGIONS_PALETTE = ct.color_palette("regions")  # Regions palette


# Module constants
_REGIONS = ["nmd_target", "distal_nmd", "start_proximal", "long_exon"]
_REGION_ORDER = ["nmd_target", "long_exon", "distal_nmd", "start_proximal"]
_REGION_LABELS = ["NMD target", "Long exon", "Distal", "Start proximal"]


# Functions
def categorise_regions(df):
    df["region"] = pd.Categorical(df["region"], categories=_REGION_ORDER, ordered=True)
    df = df.sort_values("region")
    return df


def plot_clinvar_bars(
    df, ax, height, ylabel, bar_label=True, axhline=False, xticks=False, **kwargs
):
    n = len(df)  # Number of bars
    x = np.arange(n)  # X ticks
    h = df[height]  # Height of bars (y-axis value)

    b = ax.bar(x=x / n, height=h, width=1 / (n + 1), **kwargs)

    if bar_label:
        ax.bar_label(b, fmt="%.2f")

    if xticks:
        ax.set_xticks(
            ticks=x / n,
            labels=_REGION_LABELS,
            rotation=45,
            rotation_mode="anchor",
            ha="right",
        )
    else:
        ax.set_xticks([])

    ax.tick_params(axis="x", length=0)  # Remove x ticks

    ax.set_ylabel(ylabel)

    # Optional horizontal dashed grey line
    if axhline:
        ax.axhline(y=axhline, linestyle="--", color=VIBRANT_PALETTE.grey)

    return None


# CDS footprints of NMD regions
footprint = pd.read_csv("../outputs/stats_nmd_footprint.tsv", sep="\t").pipe(
    categorise_regions
)

fig, ax = plt.subplots()
plot_clinvar_bars(footprint, ax, "footprint", "Proportion of CDS", xticks=True)

# Ascertainment in ClinVar
asrtn = pd.read_csv("../outputs/stats_clinvar_ascertainment.tsv", sep="\t").pipe(
    categorise_regions
)

fig, ax = plt.subplots()
plot_clinvar_bars(
    asrtn,
    ax,
    "prop_norm",
    "Truncating variants\nin ClinVar (normalised)",
    bar_label=False,
    axhline=1,
    xticks=True,
)

# Proportion of truncating variants in ClinVar which are VUS
clin_vus = (
    pd.read_csv("../outputs/stats_clinvar_acmg_by_region.tsv", sep="\t")
    .query("acmg == 'VUS'")
    .pipe(categorise_regions)
)

fig, ax = plt.subplots()
plot_clinvar_bars(
    clin_vus, ax, "proportion", "Proportion VUS\n(truncating variants)", xticks=True
)
