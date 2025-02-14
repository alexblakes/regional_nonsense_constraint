"""Plot ACMG classification by consequence and constraint."""

import matplotlib.pyplot as plt
import seaborn as sns

import src
from src import constants as C

FILE_IN = "data/statistics/gwas_ptvs_in_constrained_regions.tsv"
PNG = "data/plots/gwas/gwas_ptvs_in_constrained_regions.png"
SVG = "data/plots/gwas/gwas_ptvs_in_constrained_regions.svg"
logger = src.logger


def main():
    """Run as script."""

    df = src.read_data(FILE_IN)

    # Set plotting defaults
    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])
    sns.set_palette([sns.color_palette()[n] for n in [1, 0, -1]])

    # Create the plot
    fig, ax = plt.subplots(1, 1, figsize=(4 * C.CM, 6 * C.CM), layout="constrained")
    bars = ax.bar(df["constraint"], df["count"], color=sns.color_palette())
    ax.bar_label(bars)
    ax.set(ylabel="Unique PTVs")
    ax.set_xticklabels(
        labels=ax.get_xticklabels(), ha="right", rotation_mode="anchor", rotation=45
    )

    # Save plots
    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    return df


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
