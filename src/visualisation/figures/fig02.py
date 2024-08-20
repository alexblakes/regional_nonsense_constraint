"""Plot figure 2."""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src import constants as C
import src.visualisation as vis
from src.visualisation import 

PNG = "data/plots/figures/fig_02.png"
SVG = "data/plots/figures/fig_02.svg"

def main():
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_MAPS)
    fig = plt.figure(figsize=(18 * C.CM, 13 * C.CM), layout=("constrained"))
    subfigs = fig.subfigures(3, 1, height_ratios=[3, 4, 3])
    axs_top = subfigs[0].subplots(1, 2).flatten()
    ax_middle = subfigs[1].subplots(1, 1)
    axs_bottom = subfigs[2].subplots(1, 2).flatten()

    # MAPS plot
    maps = 

    # LOEUF vs OE95 plot

    # Venn

    # phyloP

    # AlphaMissense

    # Save figure
    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    pass

if __name__ == "__main__":
    main()
