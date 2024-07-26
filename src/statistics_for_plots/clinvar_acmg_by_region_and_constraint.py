"""Get ACMG classifications of PTVs in ClinVar by region and constraint."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C
from src import statistics_for_plots as sp

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/interim/clinvar_variants_constraint.tsv"
_FILE_OUT = "data/statistics/clinvar_acmg_by_region_and_constraint.tsv"
_ACMG_DICT = {"P": "P/LP", "LP": "P/LP", "VUS": "VUS", "B": "B/LB", "LB": "B/LB"}

logger = logging.getLogger(__name__)


def read_clinvar_data(path):
    return pd.read_csv(
        path,
        sep="\t",
        usecols=["csq", "region", "acmg", "constraint"],
        dtype="category",
    )


def filter_ptvs(df):
    return df.query("csq == 'stop_gained' | csq== 'frameshift_variant'").drop(
        "csq", axis=1
    )


def tidy_constraint_annotation(df):
    return (
        df.query("constraint != 'indeterminate'")
        .dropna(subset="constraint")
        .assign(constraint=lambda x: x["constraint"].cat.remove_unused_categories())
        .assign(constraint=lambda x: x["constraint"].str.capitalize())
    )


def get_acmg_proportions(df):
    return (
        df.groupby(["region", "constraint"])["acmg"]
        .value_counts(normalize=True, dropna=False)
        .rename("proportion")
    )


def write_out(df, path):
    df.to_csv(path, sep="\t")
    return df


def main():
    """Run as script."""

    clinvar = (
        read_clinvar_data(_FILE_IN)
        .pipe(filter_ptvs)
        .pipe(tidy_constraint_annotation)
        .replace({"acmg": _ACMG_DICT})
        .pipe(sp.sort_column, "acmg", ["P/LP", "VUS", "B/LB"])
        .pipe(sp.sort_column, "region")
        .sort_values(["region", "acmg"])
        .pipe(get_acmg_proportions)
        .pipe(write_out, _FILE_OUT)
    )

    return clinvar


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
