"""Plot Venn diagrams of newly constrained transcripts."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from statsmodels.stats import proportion as prop

import src
from src import constants as C
from src.constraint import plot_venn_diagrams as venn
from src.visualisation import phylop_plots as pp

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN_GNOMAD_CST = "data/final/gene_list_gnomad_constrained.txt"
_FILE_IN_NMD_TARGET = "data/final/gene_list_nmd_target_constrained.txt"
_FILE_IN_START_PROX = "data/final/gene_list_start_proximal_constrained.txt"
_FILE_IN_LONG_EXON = "data/final/gene_list_long_exon_constrained.txt"
_FILE_IN_DISTAL = "data/final/gene_list_distal_constrained.txt"
_LABELS = ["Any region", "NMD target", "Start proximal", "Long exon", "Distal"]
_SVG = "data/plots/constraint/newly_constrained.svg"
_PNG = "data/plots/constraint/newly_constrained.png"

logger = logging.getLogger(__name__)


def get_pct_new(region_set, transcript_set, label):
    total = len(region_set)
    new = len(region_set - transcript_set)

    return pd.DataFrame(
        data={
            "total_constrained": total,
            "new_constrained": new,
            "percent": (new / total) * 100,
        },
        index=[label],
    )


def get_err_95(df):
    ci_hi = prop.proportion_confint(
        df["new_constrained"],
        df["total_constrained"],
        alpha=0.05,
        method="normal",
    )[1]

    # Express the confidence interval as an error margin
    err = (ci_hi * 100) - df["percent"]

    return df.assign(err_95=err)


def main():
    """Run as script."""

    # Read gene lists
    gnomad = venn.get_gene_set(_FILE_IN_GNOMAD_CST)
    target = venn.get_gene_set(_FILE_IN_NMD_TARGET)
    start = venn.get_gene_set(_FILE_IN_START_PROX)
    long_exon = venn.get_gene_set(_FILE_IN_LONG_EXON)
    distal = venn.get_gene_set(_FILE_IN_DISTAL)
    all_regions = target | start | long_exon | distal

    # Create a list of gene sets
    regions = [all_regions, target, start, long_exon, distal]

    # Combine statistics into one dataframe
    dfs = [get_pct_new(r, gnomad, l) for r, l in zip(regions, _LABELS)]
    df = pd.concat(dfs).pipe(get_err_95).iloc[::-1] # Reverse order for plotting
    logger.info(f"Summary statistics:\n{df}")

    # Set plot style
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)
    colors = sns.color_palette()[::-1] # Reversed for plotting

    # Instantiate the figure
    fig, ax = plt.subplots(1, 1, figsize=(8.9 * C.CM, 4 * C.CM), layout="constrained")

    # Plot horizontal bars
    bars = pp.horizontal_bars(
        df["percent"],
        xerr=df["err_95"],
        color=colors,
    )
    ax.bar_label(bars, fmt="%.2f", padding=2)
    ax.set_xlabel(
        "Constrained transcripts with weak pLI\nand LOEUF scores in gnomAD v4.1 (%)"
    )

    # Save the figure
    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return fig, ax


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
