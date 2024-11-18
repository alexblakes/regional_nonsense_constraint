"""Plot Venn diagrams of newly constrained transcripts."""

import logging

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib_venn as mv
import seaborn as sns

import src
from src import constants as C

logger = logging.getLogger(__name__)


def read_data(path):
    return pd.read_csv(path, sep="\t", index_col="ensg")


def get_constrained_transcript_set(df, column):
    enst_set = set(df[df[column]].index)

    logger.info(f"Constrained transcripts ({column}): {len(enst_set)}")

    return enst_set


def main():
    """Run as script."""

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("file_in")
    parser.add_argument("target_metric", choices=["gnomad","shet"], help="Metric to compare against.")
    parser.add_argument("file_out_name")
    args = parser.parse_args()

    # Get data for Venn diagrams
    df = read_data(args.file_in)

    any_region = get_constrained_transcript_set(df, "any_region_constrained")
    nmd_target = get_constrained_transcript_set(df, "nmd_target")
    start_proximal = get_constrained_transcript_set(df, "start_proximal")
    long_exon = get_constrained_transcript_set(df, "long_exon")
    distal = get_constrained_transcript_set(df, "distal_nmd")
    gnomad = get_constrained_transcript_set(df, "gnomad_constrained")
    shet = get_constrained_transcript_set(df, "shet_constrained")

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])

    # Instantiate the figure
    fig, axs = plt.subplots(2, 3, figsize=(14 * C.CM, 7 * C.CM), layout="constrained")

    # Get parameters for Venn diagrams
    axs = list(axs.flatten())
    axs.pop(2).set_axis_off()  # The top-right Axes is left blank
    gene_sets = [any_region, nmd_target, start_proximal, long_exon, distal]
    labels = ["Any region", "NMD target", "Start proximal", "Long exon", "Distal"]
    colors = sns.color_palette()

    if args.target_metric == "gnomad":
        target_set = gnomad
        venn_label="pLI > 0.9 or\nLOEUF < 0.6"
    elif args.target_metric == "shet":
        target_set = shet
        venn_label="shet > 0.06"

    # Plot Venn diagrams
    for ax, _set, label, color in zip(axs, gene_sets, labels, colors):
        
        v = mv.venn2_unweighted(
            ax=ax,
            subsets=[target_set, _set],
            set_labels=[venn_label, label],
            set_colors=[sns.color_palette()[0], color],
        )

        v.get_label_by_id("A").set(
            color=sns.color_palette()[0], fontsize=8, ma="center"
        )
        v.get_label_by_id("B").set(color=color, fontsize=8)

    plt.savefig(f"data/plots/shet/{args.file_out_name}.png", dpi=600)
    plt.savefig(f"data/plots/shet/{args.file_out_name}.svg")
    
    plt.close("all")

    return None


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
