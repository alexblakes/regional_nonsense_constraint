"""Plot ACMG classification of nonsense / frameshift variants in ClinVar."""

import itertools
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src.clinvar.annotate_with_constraint import read_clinvar

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/interim/clinvar_variants_constraint.tsv"
_ACMG_DICT = {"P": "P/LP", "LP": "P/LP", "B": "B/LB", "LB": "B/LB"}
_FIG_OUT = "data/plots/clinvar/clinvar_acmg"

logger = logging.getLogger(__name__)


def read_clinvar_data(path):
    return (
        pd.read_csv(
            path,
            sep="\t",
            usecols=["csq", "region", "acmg", "constraint"],
        )
        .query("csq == 'stop_gained' | csq== 'frameshift_variant'")
        .query("constraint != 'indeterminate'")
        .drop("csq", axis=1)
    )


def tidy_labels(df):
    df = df.replace(
        {
            "acmg": _ACMG_DICT,
            "region": C.NMD_REGIONS_DICT,
        }
    )
    df["constraint"] = df["constraint"].str.capitalize()

    # Sort by region and ACMG
    df["region"] = pd.Categorical(
        df["region"], categories=C.NMD_REGION_LABELS, ordered=True
    )
    df["acmg"] = pd.Categorical(
        df["acmg"], categories=["P/LP", "VUS", "B/LB"], ordered=True
    )

    df.sort_values(["region", "acmg"])

    return df


def get_acmg_proportions(df):
    grouped = df.groupby(["region", "constraint"])
    proportions = grouped["acmg"].value_counts(normalize=True, dropna=False)

    # Reindex to keep categories with zero variants
    n_index_levels = proportions.index.nlevels
    index_values = [proportions.index.levels[n] for n in range(n_index_levels)]
    new_index = itertools.product(*index_values)
    proportions = proportions.reindex(new_index, fill_value=0)

    logger.info(f"ACMG proportions:\n{proportions}")

    return proportions


def plot_bars(data, ax=None, title=None, legend=0):
    if not ax:
        ax = plt.gca()

    acmg_categories = data.index.get_level_values("acmg").unique()

    for i, acmg in enumerate(acmg_categories):
        multiplier = i
        _data = data.xs(acmg, level="acmg")
        n_groups = len(_data)
        x = np.arange(n_groups)
        n_bars = len(acmg_categories)
        bar_width = 1 / (n_bars + 1)
        offset = bar_width * multiplier

        bars = ax.bar(x + offset, _data, bar_width)
        ax.set_xticks(x + bar_width, labels=_data.index)
        ax.set_ylabel("Proportion of variants")
        ax.label_outer()
        ax.set_title(title)

    if legend:
        bars = [b for b in ax.containers]
        ax.legend(bars, acmg_categories)

    return None


def main():
    """Run as script."""

    # Parse and tidy the data
    clinvar = read_clinvar_data(_FILE_IN).pipe(tidy_labels).pipe(get_acmg_proportions)

    # Style and palette
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)
    sns.set_palette([sns.color_palette()[n] for n in [1, 0, -1]])

    # Instantiate the figure
    fig, axs = plt.subplots(
        1,
        4,
        figsize=(18 * C.CM, 4 * C.CM),
        sharey=True,
        layout="constrained",
    )

    # Create plots, stratified by region
    regions = clinvar.index.get_level_values("region").unique()
    legends = [0, 0, 0, 1]
    for ax, region, legend in zip(axs, regions, legends):
        data = clinvar.xs(region, level="region")
        plot_bars(data, ax, title=region, legend=legend)

    # Save the figure
    plt.savefig(f"{_FIG_OUT}.png", dpi=600)
    plt.savefig(f"{_FIG_OUT}.svg")

    return None


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
