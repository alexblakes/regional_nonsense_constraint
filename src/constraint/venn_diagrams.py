"""Plot Venn diagrams of newly constrained transcripts."""

import logging
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib_venn as mv
import seaborn as sns

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

logger = logging.getLogger(__name__)


def main():
    """Run as script."""

    # Read gene lists
    gene_set = lambda x: set(pd.read_csv(x, header=None).iloc[:, 0].to_list())

    gnomad = gene_set(C.GENE_LIST_GNOMAD_CST)
    target = gene_set(C.GENE_LIST_NMD_TARGET)
    start = gene_set(C.GENE_LIST_START_PROX)
    long_exon = gene_set(C.GENE_LIST_LONG_EXON)
    distal = gene_set(C.GENE_LIST_DISTAL)

    # Set plot style
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)

    # Instantiate the figure
    fig, axs = plt.subplots(
        2, 2, figsize=(8.9 * C.CM, 8.9 * C.CM), layout="constrained"
    )

    # Get parameters for Venn diagrams
    axs = axs.flatten()
    gene_sets = [target, start, long_exon, distal]
    colors = sns.color_palette()[1:]

    # Plot Venn diagrams
    for ax, _set, label, color in zip(axs, gene_sets, C.NMD_REGION_LABELS, colors):
        v = mv.venn2_unweighted(
            ax=ax,
            subsets=[gnomad, _set],
            set_labels=["gnomAD", label],
            set_colors=[sns.color_palette()[0], color],
        )

        v.get_label_by_id("A").set(color=sns.color_palette()[0])
        v.get_label_by_id("B").set(color=color)

    plt.savefig("data/plots/constraint/venn.png", dpi=600)
    plt.savefig("data/plots/constraint/venn.svg")
    plt.close("all")

    return None


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
