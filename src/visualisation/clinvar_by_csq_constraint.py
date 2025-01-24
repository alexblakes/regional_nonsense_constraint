"""Plot ACMG classification by consequence and constraint."""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import statistics_for_plots as sp
from src import visualisation as vis

FILE_IN = "data/statistics/clinvar_by_csq_constraint.tsv"
PNG = "data/plots/clinvar/clinvar_acmg_by_csq_and_constraint.png"
SVG = "data/plots/clinvar/clinvar_acmg_by_csq_and_constraint.svg"

logger = src.logger


def read_data(path=FILE_IN):
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
            "".join([x, "\n", "N=", z]) for x, z in zip(old_xticklabels, x_tick_append)
        ]

        ax.set_xticks(ax.get_xticks(), new_xticklabels)

    ax.set(*args, **kwargs)

    if legend:
        ax.legend(loc="upper left")

    return ax


def main():
    """Run as script."""

    # Load data
    clinvar = read_data(FILE_IN)

    # Set plotting defaults
    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])
    sns.set_palette([sns.color_palette()[n] for n in [1, 0, -1]])

    fig, axs = plt.subplots(
        2, 2, figsize=(12 * C.CM, 10 * C.CM), layout="constrained", sharey=True
    )

    axs = axs.flatten()
    csqs = C.CSQ_LABELS.values()
    legends = [1, 1, 1, 1]

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
            data.xs("Constrained", level="constraint")["total"].astype(str).to_list()
        )
        customise_plot(
            ax,
            legend,
            x_tick_append=variant_counts,
            title=csq,
            ylabel="Proportion of variants in ClinVar",
        )

    plt.savefig(PNG, dpi=600)
    plt.savefig(SVG)
    plt.close("all")

    return clinvar


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
