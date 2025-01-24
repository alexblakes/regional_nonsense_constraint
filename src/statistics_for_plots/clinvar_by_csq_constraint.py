"""Get constraint annotation for variants split by consequence and ACMG class.

Usage:
    python3 -m src.statistics_for_plots.clinvar_by_csq_constraint <file_in> <file_out>
"""

import argparse

import pandas as pd
from statsmodels.stats import proportion

import src
from src import statistics_for_plots as sp
from src import constants as C

# FILE_IN = "data/interim/clinvar_variants_constraint.tsv.gz"
# FILE_OUT = "data/statistics/clinvar_by_csq_constraint.tsv"
ACMG_ORDER = ["P/LP", "VUS", "B/LB"]

logger = src.logger


@src.log_step
def read_data(path):
    keep_cols = "csq acmg constraint region".split()
    return pd.read_csv(path, sep="\t", usecols=keep_cols)


def get_statistics(df):
    grouping = ["csq", "acmg"]
    return (
        df.groupby(grouping)["constraint"]
        .value_counts()
        .to_frame("count")
        .assign(
            total=lambda x: x.groupby(grouping)["count"].transform("sum"),
            proportion=lambda x: x["count"] / x.groupby(grouping)["count"].sum(),
            ci95lo=lambda x: proportion.proportion_confint(x["count"], x["total"])[0],
            err95=lambda x: x.proportion - x.ci95lo,
        )
        .drop("ci95lo", axis=1)
    )


def sort_labels(df):
    return (
        df.reset_index()
        .pipe(sp.sort_column, "csq", C.CSQ_LABELS)
        .pipe(sp.sort_column, "acmg", ACMG_ORDER)
        .pipe(sp.sort_column, "constraint", C.CONSTRAINT_LABELS)
        .sort_values(["csq", "acmg", "constraint"])
    )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("file_in", help="Input file")
    parser.add_argument("file_out", help="Output file")

    return parser.parse_args()


def main():
    """Run as script."""
    args = parse_args()

    return (
        read_data(args.file_in)
        .pipe(get_statistics)
        .pipe(sort_labels)
        .pipe(src.write_out, args.file_out)
    )


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
