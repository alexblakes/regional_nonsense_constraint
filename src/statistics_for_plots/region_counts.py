"""Count the number of transcripts with each NMD region."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import statistics_for_plots as sp

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/final/nmd_annotations_simple.tsv.gz"
_FILE_OUT = "data/statistics/region_counts.tsv"

logger = logging.getLogger(__name__)


def read_data(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "pos", "enst", "region"],
        usecols=["enst", "region"],
        dtype="category",
    )


def count_total_transcripts(df):
    return pd.Series(
        data=df.enst.nunique(), index=["transcript"], name="transcript_count"
    )


def count_transcripts_per_region(df):
    return df.groupby("region")["enst"].nunique().rename("transcript_count")


def combine_totals(df):
    return pd.concat([count_total_transcripts(df), count_transcripts_per_region(df)])


def write_out(series, path):
    series.to_csv(path, sep="\t")
    return series


def main():
    """Run as script."""

    df = (
        read_data(_FILE_IN)
        .pipe(combine_totals)
        .pipe(sp.sort_index)
        .pipe(write_out, _FILE_OUT)
    )

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
