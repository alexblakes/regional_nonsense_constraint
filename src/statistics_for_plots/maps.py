"""Get plot statistics for MAPS."""

import logging

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import src
from src import constants as C
from src import statistics_for_plots as sp

_FILE_IN = "data/interim/maps.tsv"
_FILE_OUT = "data/statistics/maps.tsv"

logger = logging.getLogger(__name__)


def read_data(path=_FILE_IN):
    return pd.read_csv(path, sep="\t")


def select_rows(df):
    m1 = df["csq"] == "stop_gained"
    m2 = df["region"] == "transcript"

    return df[m1 | m2]


def tidy_maps(df):
    df = sp.sort_column(df, "csq", C.CSQS)
    df = sp.sort_column(df, "region", C.REGIONS)

    labels = pd.Categorical(C.MAPS_LABELS, C.MAPS_LABELS, ordered=True)
    df = (
        df.sort_values(["csq", "region"])
        .assign(csq=labels)
        .sort_values("csq", ascending=False)
    )

    return df


def write_out(df, path):
    df.to_csv(path, sep="\t")
    return df


def main():
    """Run as script."""

    return read_data().pipe(select_rows).pipe(tidy_maps).pipe(write_out, _FILE_OUT)


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
