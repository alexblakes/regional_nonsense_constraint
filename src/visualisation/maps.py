"""Create a horizontal point-range plot for MAPS scores."""

import logging

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import src
from src import constants as C

_FILE_IN = "data/statistics/maps.tsv"
_PNG = "data/plots/maps/maps.png"
_SVG = "data/plots/maps/maps.svg"

logger = logging.getLogger(__name__)


def read_data(path=_FILE_IN):
    return pd.read_csv(path, sep="\t", index_col="csq")


def h_pointrange(y, x, xerr, ax=None):
    if ax == None:
        ax = plt.gca()

    c = sns.color_palette()

    ax.scatter(y=y, x=x, c=c)
    ax.errorbar(y=y, x=x, xerr=xerr, ecolor=c, linestyle="")

    return None

def customise_plot(ax, **kwargs):
    if not ax:
        ax = plt.gca()

    ax.set_xlabel("MAPS")
    
    ax.set(**kwargs)

    return ax    


def plot(ax=None):
    if not ax:
        ax = plt.gca()

    maps = read_data()
    h_pointrange(y=maps.index, x=maps["maps"], xerr=maps["ci95"], ax=ax)
    customise_plot(ax)

    return ax


def main():
    """Run as script."""

    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_MAPS)
    fig, ax = plt.subplots(1, 1, figsize=(8 * C.CM, 4 * C.CM), layout="constrained")

    plot()

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return maps


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
