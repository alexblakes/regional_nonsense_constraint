"""Plot extended data figure 5"""

import logging

from matplotlib import pyplot as plt

import src
from src import constants as C
from src import visualisation as vis
from src.visualisation import orthogonal_metrics_bxp as omb

PNG = "data/plots/figures/ext_fig_05.png"
SVG = "data/plots/figures/ext_fig_05.svg"

logger = logging.getLogger(__name__)


def main():
    """Run as script."""

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])
    fig, axs = plt.subplots(2, 2, figsize=(12 * C.CM, 12 * C.CM), layout="constrained")
    axs = axs.flatten()

    omb.plot(axs)

    # Add panel labels
    labels = list("abcd")
    for ax, label in zip(axs, labels):
        vis.panel_label(ax, label)

    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    return fig, axs


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
