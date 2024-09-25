"""Plot figure 1."""

import pandas as pd
import matplotlib.pyplot as plt

from src import constants as C
from src import visualisation as vis
from src.visualisation import regions_cds_proportions as rcp
from src.visualisation import clinvar_ptv_ascertainment as cpa
from src.visualisation import clinvar_proportion_vus_by_region as cpv

CDS_PROPORTION = "data/statistics/regions_cds_proportions.tsv"
CLINVAR_ASCERTAINMENT = "data/statistics/clinvar_ptv_ascertainment.tsv"
CLINVAR_VUS = "data/statistics/clinvar_vus_by_region.tsv"
NMD_DIAGRAM = "data/plots/manual/fig_1_transcript_diagram.png"
G_ABSTRACT = "data/plots/manual/fig_1_graphical_abstract.png"
PNG = "data/plots/figures/fig_01.png"
SVG = "data/plots/figures/fig_01.svg"


def main():
    # Get data
    graphical_abstract = plt.imread(G_ABSTRACT)
    transcript_diagram = plt.imread(NMD_DIAGRAM)

    # Plotting style
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)

    # Instantiate the figure
    fig = plt.figure(
        figsize=(12 * C.CM, 12 * C.CM),
        layout="constrained",
    )
    subfigs = fig.subfigures(1, 2, width_ratios=(7, 5))
    ax_left_top, ax_left_bottom = (
        subfigs[0].subplots(2, 1, height_ratios=(1, 2), gridspec_kw=dict(hspace=0.1)).flatten()
    )
    axs_right = subfigs[1].subplots(3, 1).flatten()

    # Biorender graphical abstract
    ax_left_top.imshow(graphical_abstract)
    ax_left_top.axis("off")

    # Biorender transcript diagram
    ax_left_bottom.imshow(transcript_diagram)
    ax_left_bottom.axis("off")

    # NMD region footprints
    cds_proportion = rcp.read_cds_proportions(CDS_PROPORTION)
    vis.vertical_bars(cds_proportion, ax=axs_right[0])
    rcp.customise_plot(ax=axs_right[0])

    # Ascertainment
    ptv_ascertainment = cpa.read_ptv_ascertainment(CLINVAR_ASCERTAINMENT)
    vis.vertical_bars(ptv_ascertainment, ax=axs_right[1])
    cpa.customise_plot(ax=axs_right[1])

    # Proportion VUS
    proportion_vus = cpv.read_vus_proportions(CLINVAR_VUS)
    vis.vertical_bars(
        proportion_vus["proportion_vus"], ax=axs_right[2], yerr=proportion_vus["err"]
    )
    cpv.customise_plot(proportion_vus, ax=axs_right[2])

    # Tidy axes
    for ax in axs_right:
        ax.label_outer()

    # Panel labels
    axes = [ax_left_top, ax_left_bottom] + list(axs_right)
    labels = list("abcde")
    xs = [0.05] * 2 + [-0.05] * 3
    ys = [0.98] *2 + [1.05] * 3

    for ax, label, x, y in zip(axes, labels, xs, ys):
        vis.panel_label(ax, label, x, y)

    # Save the figure
    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    pass


if __name__ == "__main__":
    main()
