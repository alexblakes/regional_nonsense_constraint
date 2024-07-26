"""Plot ACMG classification of nonsense / frameshift variants in ClinVar."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/statistics/clinvar_acmg_by_region_and_constraint.tsv"
_PNG = "data/plots/clinvar/clinvar_acmg_by_region_and_constraint.png"
_SVG = "data/plots/clinvar/clinvar_acmg_by_region_and_constraint.svg"
_REGION_LABELS = ["NMD target", "Start proximal", "Long exon", "Distal"]

logger = logging.getLogger(__name__)


def read_data(path):
    return pd.read_csv(
        path, sep="\t", index_col=["region", "constraint", "acmg"]
    ).squeeze()


def customise_plot(ax=None, title=None, legend=False):
    if not ax:
        ax = plt.gca()

    ax.tick_params(labelbottom=True)
    ax.set_ylabel("Proportion of PTVs in ClinVar")
    ax.label_outer()

    ax.set_title(title)

    if legend:
        ax.legend()

    return ax


def main():
    """Run as script."""

    clinvar = read_data(_FILE_IN)

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])
    sns.set_palette([sns.color_palette()[n] for n in [1, 0, -1]])

    fig, axs = plt.subplots(
        1, 4, figsize=(18 * C.CM, 4 * C.CM), sharey=True, layout="constrained"
    )

    axs = axs.flatten()
    regions = _REGION_LABELS
    legends = [0, 0, 0, 1]

    for ax, region, legend in zip(axs, regions, legends):
        data = clinvar.xs(region, level="region")
        vis.vertical_grouped_bars(data, ax)
        customise_plot(ax, region, legend)

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    # plt.close("all")

    return axs


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
