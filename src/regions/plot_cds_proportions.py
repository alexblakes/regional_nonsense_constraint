"""Plot ACMG classification of nonsense / frameshift variants in ClinVar."""

from locale import normalize
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/final/nmd_annotations_simple.tsv.gz"
_FILE_OUT = "data/statistics/regions_cds_proportions.tsv"
_PNG = "data/plots/regions_cds_proportions.png"
_SVG = "data/plots/regions_cds_proportions.svg"
_CATEGORIES = ["full_cds", "nmd_target", "start_proximal", "long_exon", "distal_nmd"]
_LABELS = ["Full CDS", "NMD target", "Start proximal", "Long exon", "Distal"]

logger = logging.getLogger(__name__)


def read_nmd_regions(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "pos", "enst", "region"],
        nrows=100000,
        usecols=["region"],
        dtype="category",
    ).loc[:, "region"]


def proportion_value_counts(series):
    return series.value_counts(normalize=True).rename("proportion")


def append_full_cds(series):
    """Append the proportion of sites in the full CDS (1)."""
    cds = pd.Series(data=[1.0], index=["full_cds"], name="proportion")
    return pd.concat([series, cds])


def sort_index(series, categories, labels, **kwargs):
    """Create an ordered categorical index for a series."""
    _ix = pd.CategoricalIndex(categories, **kwargs)
    return series.reindex(_ix).rename(dict(zip(categories, labels)))


def write_out(series, path):
    series.to_csv(path, sep="\t")
    return series


def plot_cds_proportions(series, ax=None):
    if not ax:
        ax = plt.gca()

    colors = sns.color_palette()
    x = series.index

    bars = ax.bar(
        x,
        height=series,
        color=colors,
    )

    ax.bar_label(bars, fmt="{:.2%}")
    ax.set_xticks(ticks=x, labels=x, rotation=45, ha="right", rotation_mode="anchor")

    return ax


def main():
    """Run as script."""
    region_counts = (
        read_nmd_regions(_FILE_IN)
        .pipe(proportion_value_counts)
        .pipe(append_full_cds)
        .pipe(sort_index, _CATEGORIES, _LABELS, ordered=True, name="region")
        .pipe(write_out, _FILE_OUT)
    )

    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)

    fig, ax = plt.subplots(1, 1, figsize=(8.9 * C.CM, 8.9 * C.CM), layout="constrained")

    plot_cds_proportions(region_counts)

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)

    return region_counts


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
