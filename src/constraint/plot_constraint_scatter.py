"""Plot constraint pairplots."""

import itertools
import logging
from pathlib import Path

import adjustText
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import seaborn as sns

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_REGIONAL_NONSENSE_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
_GENE_IDS = "data/interim/gene_ids.tsv"
_PNG = "data/plots/constraint/region_pair_plots.png"
_SVG = "data/plots/constraint/region_pair_plots.svg"
_REGION_NAMES = ["nmd_target", "start_proximal", "long_exon", "distal_nmd"]
_REGION_LABELS = ["NMD target", "Start proximal", "Long exon", "Distal"]

logger = logging.getLogger(__name__)


def read_constraint(path=_REGIONAL_NONSENSE_CONSTRAINT):
    return pd.read_csv(
        path,
        sep="\t",
        usecols=[
            "enst",
            "region",
            "n_exp",
            "oe",
            "oe_ci_hi",
            "p",
            "pli",
            "loeuf",
            "syn_p",
            "constraint",
        ],
    )


def filter_constraint(df):
    return (
        df.dropna()
        .loc[lambda x: x["syn_p"] >= stats.norm.cdf(-1), :]
        .loc[lambda x: x["n_exp"] >= 5, :]
    )


def get_gene_ids(path=_GENE_IDS):
    return pd.read_csv(
        path,
        sep="\t",
        usecols=["transcript_id", "gene_name"],
        header=0,
    ).set_axis(["enst", "symbol"], axis=1)


def parse_data():
    constraint = read_constraint().pipe(filter_constraint)
    gene_ids = get_gene_ids()

    return constraint.merge(gene_ids, validate="many_to_one")


def reorient_data(df):
    return (
        df.pivot(
            index=["enst", "symbol"],
            columns="region",
            values=["constraint", "oe_ci_hi"],
        )
        .reorder_levels([1, 0], axis=1)
        .sort_index(axis=1)
    )


def rename_columns(df):
    df.columns = ["_".join([x, y]) for x, y in df.columns]
    return df


def keep_relevant_columns(df, xlabel, ylabel):
    # Get column names
    x_constraint = xlabel + "_constraint"
    y_constraint = ylabel + "_constraint"
    x_oe = xlabel + "_oe_ci_hi"
    y_oe = ylabel + "_oe_ci_hi"

    return (
        df[[x_constraint, x_oe, y_constraint, y_oe]]
        .copy()
        .reset_index("symbol")
        .rename(
            columns={
                x_constraint: "x_constraint",
                y_constraint: "y_constraint",
                x_oe: "x_oe_ci_hi",
                y_oe: "y_oe_ci_hi",
            }
        )
    )


def keep_constrained_regions(df):
    m1 = df["x_constraint"] == "constrained"
    m2 = df["y_constraint"] == "constrained"

    return df.loc[m1 | m2, :]


def main():
    """Run as script."""

    df = parse_data().pipe(reorient_data).pipe(rename_columns)

    # Instantiate figure
    plt.style.use(C.STYLE_DEFAULT)
    plt.style.use(C.COLOR_REGIONS)

    fig, axs = plt.subplots(
        4,
        4,
        figsize=(18 * C.CM, 18 * C.CM),
        layout="tight",
        sharex="col",
        sharey="row",
        subplot_kw=dict(box_aspect=1),
    )

    # Define iterables for plotting
    axs = axs.ravel("F")
    subsets = itertools.product(_REGION_NAMES, _REGION_NAMES)
    palette = sns.color_palette()
    x_colors = [palette[i] for i in [1] * 4 + [2] * 4 + [3] * 4 + [4] * 4]
    y_colors = [palette[i] for i in [1, 2, 3, 4] * 4]
    grey = palette[0]
    region_dict = {a: b for a, b in zip(_REGION_NAMES, _REGION_LABELS)}

    # Create plots
    for ax, (xlabel, ylabel), x_color, y_color in zip(axs, subsets, x_colors, y_colors):
        ax.set_xlabel(f"O/E upper 95% CI\n{region_dict[xlabel]}")
        ax.set_ylabel(f"O/E upper 95% CI\n{region_dict[ylabel]}")
        ax.label_outer()

        if xlabel == ylabel:
            continue

        # Subset the data for each scatter plot
        ax_data = keep_relevant_columns(df, xlabel, ylabel).pipe(
            keep_constrained_regions
        )

        # Plot the data in three parts
        m1 = ax_data["x_constraint"] == "constrained"
        m2 = ax_data["y_constraint"] == "constrained"

        x_only = ax_data[m1 & ~m2]
        y_only = ax_data[~m1 & m2]
        both = ax_data[m1 & m2]

        for data, color in zip(
            [x_only, y_only, both], [x_color, y_color, grey]
        ):
            # Plot scatter plots
            ax.scatter(
                data["x_oe_ci_hi"],
                data["y_oe_ci_hi"],
                alpha=0.5,
                linewidth=0,
                color=color,
            )

        # Highlight the strongest outliers
        ax_data["diff"] = (ax_data["x_oe_ci_hi"] - ax_data["y_oe_ci_hi"]).astype(float)
        min3 = ax_data.nsmallest(5, "diff")

        ax.scatter(
            min3["x_oe_ci_hi"],
            min3["y_oe_ci_hi"],
            color="None",
            edgecolor="black",
            linewidth=0.5,
        )

        # Annotate the N largest outliers
        annots = []
        for i, row in min3.iterrows():
            annots.append(
                ax.text(
                    x=row["x_oe_ci_hi"],
                    y=row["y_oe_ci_hi"],
                    s=row["symbol"],
                    ha="left",
                    va="top",
                    size=7,
                )
            )
        adjustText.adjust_text(
            annots,
            ax=ax,
            expand=(1.1, 1.2),
            time_lim=3,
            expand_axes=False,
            only_move="x+y",
            avoid_self=True,
            arrowprops=dict(arrowstyle="-", linewidth=0.5),
        )

    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
