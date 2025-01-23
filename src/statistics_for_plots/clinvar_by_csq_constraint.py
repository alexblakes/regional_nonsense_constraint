"""Get constraint annotation for variants split by consequence and ACMG class."""

import pandas as pd
from statsmodels.stats import proportion

import src
from src import statistics_for_plots as sp
from src import constants as C

FILE_IN = "data/interim/clinvar_variants_constraint.tsv.gz"
FILE_OUT = "data/statistics/clinvar_by_csq_constraint.tsv"
ACMG_ORDER = ["P/LP", "VUS", "B/LB"]

logger = src.logger


@src.log_step
def read_data(path=FILE_IN):
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
        .sort_values(["csq","acmg","constraint"])
    )


def main():
    """Run as script."""
    return (
        read_data()
        .pipe(get_statistics)
        .pipe(sort_labels)
        .pipe(src.write_out, FILE_OUT)
    )


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
