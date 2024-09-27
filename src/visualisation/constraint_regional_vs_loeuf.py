"""Create a hexbin plot for regional constraint vs LOEUF scores."""

import logging
from re import I

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

import src
from src import constants as C

FILE_IN = "data/final/regional_nonsense_constraint.tsv"
PNG = "data/plots/constraint/oe_ci_hi_vs_loeuf.png"
SVG = "data/plots/constraint/oe_ci_hi_vs_loeuf.svg"
logger = logging.getLogger(__name__)


def read_data(path=FILE_IN):
    return pd.read_csv(path, sep="\t")

def filter_transcripts(df):
    return df.query("region == 'transcript'")

def hexbin(x, y, ax=None, colorbar=True, **kwargs):
    kwargs.setdefault("mincnt", 1)
    kwargs.setdefault("gridsize", 40)
    kwargs.setdefault("linewidth", 0)
    kwargs.setdefault("cmap", "Greys")

    if not ax:
        ax = plt.gca()

    im = ax.hexbin(x, y, **kwargs)

    if colorbar:
        plt.colorbar(im, ax=ax, label="Transcripts", shrink=0.8, aspect=15)

    return None


def customise_plot(rho, ax=None):
    if not ax:
        ax = plt.gca()

    # Axis labels
    ax.set_xlabel("Nonsense OE95")
    ax.set_ylabel("LOEUF")

    # Add rho annotation
    ax.text(
        x=0.95,
        y=0.05,
        s=f"rho={rho:.3}",
        ha="right",
        va="bottom",
        transform=ax.transAxes,
    )

    # Aspect
    ax.set_box_aspect(1)

    return ax

def plot(ax=None):
    if not ax:
        ax = plt.gca()

    df = read_data().pipe(filter_transcripts)
    rho, pvalue = stats.spearmanr(df.oe_ci_hi, df.loeuf, nan_policy="omit") # type: ignore
    logger.info(f"Spearman rho: {rho}, pvalue: {pvalue}")
    hexbin(df.oe_ci_hi, df.loeuf, extent=[0,2,0,2], ax=ax)
    customise_plot(rho, ax)

    return ax


def main():
    """Run as script."""

    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)
    fig, ax = plt.subplots(1, 1, figsize=(8.9 * C.CM, 8.9 * C.CM))

    plot()

    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    return fig, ax


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
