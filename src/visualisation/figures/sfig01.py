"""Plot Supplementary Figure 1."""

import logging

from matplotlib import pyplot as plt
import pandas as pd

import src
from src import constants as C
from src import visualisation as vis
from src.visualisation import constraint_oe_plots as cop
from src.visualisation import constraint_regional_vs_loeuf as rvl


PNG = "data/plots/figures/supp_fig_01.png"
SVG = "data/plots/figures/supp_fig_01.svg"

logger = logging.getLogger(__name__)


def main():
    """Run as script."""

    plt.style.use(C.STYLE_DEFAULT)

    fig = plt.figure(figsize=(18 * C.CM, 12 * C.CM), layout="constrained")
    subfig_top, subfig_bottom = fig.subfigures(2,1)
    axs_top = subfig_top.subplots(1,3)
    axs_bottom = subfig_bottom.subplots(1,3)
    axs_bottom[1].axis("off")
    axs_bottom[2].axis("off")

    cop.plot(subfig_top, axs_top)
    rvl.plot(axs_bottom[0])

    axs = list(axs_top) + [axs_bottom[0]]
    labels = list("abcd")
    for ax, label in zip(axs, labels):
        vis.panel_label(ax, label)

    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    return None


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
