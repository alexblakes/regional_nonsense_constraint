"""Module docstring."""

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

def plot_venn(ax=None):
    if not ax: 
        ax = plt.gca()
    
    mv.venn2(
        ax=ax,
        subsets=[gnomad, target],
        set_labels=["gnomAD", "NMD target"],
        set_colors=[palette[0], palette[1]]
    )

    return ax

def main():
    """Run as script."""

    # Set plot style
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)

    # Read gene lists
    gene_set = lambda x: set(pd.read_csv(x, header=None).iloc[:, 0].to_list())

    gnomad = gene_set(C.GENE_LIST_GNOMAD_CST)
    target = gene_set(C.GENE_LIST_NMD_TARGET)
    start = gene_set(C.GENE_LIST_START_PROX)
    long_exon = gene_set(C.GENE_LIST_LONG_EXON)
    distal = gene_set(C.GENE_LIST_DISTAL)

    fig, axs = plt.subplots(2, 2, figsize=(8.9 * C.CM, 8.9 * C.CM), layout="constrained")
    
    axs = axs.flatten()
    gene_sets = [target, start, long_exon, distal]
    gene_labels = ["NMD target","Start proximal","Long exon","Distal"]
    colors = sns.color_palette()[1:]

    for ax, _set, label, color in zip(axs, gene_sets, gene_labels, colors):
        v = mv.venn2_unweighted(
            ax=ax,
            subsets=[gnomad, _set],
            set_labels=["gnomAD", label],
            set_colors=[sns.color_palette()[0], color],
        )
        
        v.get_label_by_id("A").set(color=sns.color_palette()[0])
        v.get_label_by_id("B").set(color=color)

    plt.savefig("data/plots/constraint/venn.png", dpi=600)
    plt.close("all")

    return None


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
