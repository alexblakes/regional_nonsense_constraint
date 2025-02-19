"""Plot ACMG classification of nonsense / frameshift variants in ClinVar."""



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import src
from src import constants as C
from src import visualisation as vis


_FILE_IN = "data/statistics/clinvar_ptv_ascertainment.tsv"
_PNG = "data/plots/clinvar/clinvar_ptvs_ascertainment.png"
_SVG = "data/plots/clinvar/clinvar_ptvs_ascertainment.svg"

logger = src.logger


def read_ptv_ascertainment(path):
    return pd.read_csv(path, sep="\t", index_col="region").squeeze()


def customise_plot(ax=None):
    if not ax:
        ax = plt.gca()

    # x axis
    labels = ax.get_xticklabels()
    ticks = np.arange(len(labels))
    ax.set_xticks(ticks, labels, rotation=45, ha="right", rotation_mode="anchor")

    # y axis
    ax.set_ylabel("Nonsense and frameshift\nvariant counts in ClinVar\n(normalised)")

    # Bar labels
    bars = ax.containers[0]
    ax.bar_label(bars, fmt="{:.2f}", padding=5)

    # Hline
    ax.axhline(1, color="grey", ls="--")

    # Significance asterisks
    vis.add_significance_asterisk(xs=[1,3], ys=[1.18, 1.11], ps=[1,1], ax=ax, y_adj=3)
    vis.add_significance_asterisk(xs=[2,4], ys=[0.63, 0.6], ps=[1,1], ax=ax, y_adj=3, marker="$\u2020$")

    return ax


def main():
    """Run as script."""
    ptv_ascertainment = read_ptv_ascertainment(_FILE_IN)

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])

    fig, ax = plt.subplots(1, 1, figsize=(5 * C.CM, 5 * C.CM), layout="constrained")

    vis.vertical_bars(ptv_ascertainment)
    customise_plot()

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)

    return ax


if __name__ == "__main__":
    src.add_log_handlers()
    main()
