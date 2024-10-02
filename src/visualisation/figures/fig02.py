"""Plot figure 2."""

import matplotlib.pyplot as plt

from src import constants as C
import src.visualisation as vis
from src.visualisation import maps
from src.visualisation import constraint_regional_vs_loeuf as rvl
from src.visualisation import orthogonal_metrics_bxp as omb
from src.visualisation import constraint_oe95_kde as kdes

VENN = "data/plots/constraint/venn.png"
PNG = "data/plots/figures/fig_02.png"
SVG = "data/plots/figures/fig_02.svg"

def main():

    # Instantiate figure
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_MAPS)
    fig = plt.figure(figsize=(12 * C.CM, 22 * C.CM), layout=("constrained"))
    subfigs = fig.subfigures(4, 1, height_ratios=[5, 4, 8, 10], hspace=0.05) #type: ignore
    axs_top = subfigs[0].subplots(1, 2).flatten()
    axs_middle_upper = subfigs[1].subplots(1, 5).flatten()
    ax_middle_lower = subfigs[2].subplots(1, 1)
    axs_bottom = subfigs[3].subplots(2, 2, sharex="col").flatten()

    # MAPS
    maps.plot(ax=axs_top[0])
        
    # # Nonsense OE95 vs LOEUF
    # rvl.plot(ax=axs_top[1])

    # Blank top right axes
    axs_top[1].axis("off")

    # OE95 KDE plots
    plt.style.use(C.COLOR_REGIONS)
    kdes.plot(axs_middle_upper)
    plt.style.use(C.COLOR_MAPS)

    # Venn
    venn = plt.imread(VENN)
    ax_middle_lower.imshow(venn)
    ax_middle_lower.axis("off")

    # Orthogonal metrics
    plt.style.use(C.COLOR_REGIONS)
    omb.plot(axs_bottom)

    # Panel labels
    axs = [axs_top[0], axs_middle_upper[0], ax_middle_lower] + list(axs_bottom)
    labels = list("abcdefg")
    xs = [-0.05] * 2 + [0.05] + [-0.05] * 4
    ys = [1.05] * 2 + [0.95] + [1.05] * 4
    for ax, label, x, y in zip(axs, labels, xs, ys):
        vis.panel_label(ax, label, x, y) 

    # Save figure
    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    pass

if __name__ == "__main__":
    main()
