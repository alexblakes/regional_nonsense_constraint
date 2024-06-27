"""Create a horizontal point-range plot for MAPS scores."""

import logging
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import src
from src import constants as C

_FILE_IN = "data/interim/maps.tsv"
_PNG = "data/plots/maps/maps.png"
_SVG = "data/plots/maps/maps.svg"
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

logger = logging.getLogger(__name__)


def read_data(path=_FILE_IN):
    return pd.read_csv(path, sep="\t")


def select_rows(df):
    m1 = df["csq"] == "stop_gained"
    m2 = df["region"] == "transcript"

    return df[m1 | m2]


def tidy_maps(df):
    df["csq"] = pd.Categorical(df["csq"], C.CSQS, ordered=True)
    df["region"] = pd.Categorical(df["region"], C.REGIONS, ordered=True)

    _labels = pd.Categorical(C.MAPS_LABELS, C.MAPS_LABELS, ordered=True)
    df = (
        df.sort_values(["csq", "region"])
        .assign(csq=_labels)
        .sort_values("csq", ascending=False)
    )

    return df


def plot_maps(y, x, xerr, ax=None):
    # Get current axis if not specified
    if ax == None:
        ax = plt.gca()

    # Point range plot
    _c = sns.color_palette()
    ax.scatter(y=y, x=x, c=_c)
    ax.errorbar(y=y, x=x, xerr=xerr, ecolor=_c, linestyle="None")

    return None


def main():
    """Run as script."""

    df = read_data().pipe(select_rows).pipe(tidy_maps)

    # Plot style
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_MAPS)

    fig, ax = plt.subplots(1, 1, figsize=(8 * C.CM, 4 * C.CM), layout="constrained")

    plot_maps(y=df["csq"], x=df["maps"], xerr=df["ci95"])

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)


    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
