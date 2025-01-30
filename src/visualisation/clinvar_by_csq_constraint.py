"""Plot ACMG classification by consequence and constraint."""

import argparse
import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis

logger = src.logger


def read_data(path):
    index_cols = ["csq", "constraint", "acmg"]
    return pd.read_csv(path, sep="\t", index_col=index_cols)


def customise_plot(ax=None, legend=False, x_tick_append=None, *args, **kwargs):
    if not ax:
        ax = plt.gca()

    ax.tick_params(labelbottom=True, labelleft=True)

    # Add variant counts to xticklabels
    if x_tick_append:
        old_xticklabels = [x.get_text() for x in ax.get_xticklabels()]
        new_xticklabels = [
            "".join([x, "\n", "N=", z])
            for x, z in itertools.zip_longest(
                old_xticklabels, x_tick_append, fillvalue="0"
            )
        ]

        ax.set_xticks(ax.get_xticks(), new_xticklabels)

    ax.set(*args, **kwargs)

    if legend:
        ax.legend()
        # ax.legend(loc="upper left")

    return ax


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("file_in", help="Input file")
    parser.add_argument("png", help=".png output")

    return parser.parse_args()


def main():
    """Run as script."""

    args = parse_args()

    # Load data
    clinvar = read_data(args.file_in)

    # Set plotting defaults
    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])
    sns.set_palette([sns.color_palette()[n] for n in [1, 0, -1]])

    # Instantiate figure
    fig, axs = plt.subplots(
        2, 2, figsize=(12 * C.CM, 10 * C.CM), layout="constrained", sharey=True
    )

    axs = axs.flatten()
    csqs = C.CSQ_LABELS.values()
    legends = [1, 1, 1, 1]

    # Plot individual axes
    for ax, csq, legend in zip(axs, csqs, legends):
        data = clinvar.xs(csq, level="csq")
        vis.vertical_grouped_bars_with_errorbars(
            data=data,
            data_column="proportion",
            yerr_column="err95",
            bar_group="constraint",
            ax=ax,
        )

        variant_counts = (
            data.groupby(level="acmg", sort=False)["count"].sum().astype(str).to_list()
        )

        customise_plot(
            ax,
            legend,
            x_tick_append=variant_counts,
            title=csq,
            ylabel="Proportion of variants in ClinVar",
        )

    for ax, text in zip(axs, "abcd"):
        vis.panel_label(ax, text)

    # Save plots
    png = args.png
    svg = Path(png).with_suffix(".svg")
    plt.savefig(png, dpi=600)
    plt.savefig(svg)
    plt.close("all")

    return clinvar


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
