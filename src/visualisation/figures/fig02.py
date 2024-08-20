"""Plot figure 2."""

import matplotlib.pyplot as plt
import seaborn as sns

from src import constants as C
import src.visualisation as vis
from src.visualisation import maps
from src.visualisation import constraint_regional_vs_loeuf as rvl
from src.visualisation import orthogonal_metrics as om

VENN = "data/plots/constraint/venn.png"
PNG = "data/plots/figures/fig_02.png"
SVG = "data/plots/figures/fig_02.svg"

def main():

    # Instantiate figure
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_MAPS)
    fig = plt.figure(figsize=(12 * C.CM, 14 * C.CM), layout=("constrained"))
    subfigs = fig.subfigures(3, 1, height_ratios=[5, 7, 7]) #type: ignore
    axs_top = subfigs[0].subplots(1, 2).flatten()
    ax_middle = subfigs[1].subplots(1, 1)
    axs_bottom = subfigs[2].subplots(1, 2).flatten()

    # MAPS
    maps.plot(ax=axs_top[0])
    
    # Nonsense OE95 vs LOEUF
    rvl.plot(ax=axs_top[1])

    # Venn
    venn = plt.imread(VENN)
    ax_middle.imshow(venn)
    ax_middle.axis("off")

    # Orthogonal metrics
    plt.style.use(C.COLOR_REGIONS)
    om_data = om.read_data().pipe(om.sort_data)

    metrics = ["phylop", "alpha_mis"]
    labels = ["phyloP", "AlphaMissense"]
    legends = [True, False]

    for ax, metric, label, legend in zip(axs_bottom, metrics, labels, legends):
        om.grouped_boxplot(om_data, ax, y=metric)
        om.customise_plot(ax, ylabel=label, legend=legend)
        om.recolor(ax, sns.color_palette())

    # Panel labels
    axs = list(axs_top) + [ax_middle] + list(axs_bottom)
    labels = list("ABCDE")
    xs = [-0.05] * 2 + [0.05] + [-0.05] * 2
    ys = [1.05] * 2 + [0.95] + [1.05] * 2
    for ax, label, x, y in zip(axs, labels, xs, ys):
        vis.panel_label(ax, label, x, y) 

    # Save figure
    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    pass

if __name__ == "__main__":
    main()
