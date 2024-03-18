"""Summary statistics for CADD scores in constrained regions."""

# Imports
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import src
from src import constants as C
from src import statistics_for_plots


# Module constants
_DTYPE = {"pos": np.int32, "cadd_phred": np.float16}
_USECOLS = ["csq", "cadd_phred", "region", "constraint"]
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"


# Logging
logger = logging.getLogger(__name__)


# Functions
def read_cadd_annotations(path):
    logger.info("Reading annotated CADD scores.")

    return (
        pd.read_csv(
            path,
            sep="\t",
            dtype=_DTYPE,
            usecols=_USECOLS,
            low_memory=False,
            # nrows=1000000,
        )
        .dropna(subset="constraint")
        .replace({"region": "distal_nmd"}, value="distal")
    )


def get_cadd_stats(df, groupby=["constraint", "region"]):
    """Get summary stats of CADD scores by region and constraint."""

    # Lambda functions for 95% confidence intervals
    ci_l = lambda x: x["mean"] - 1.96 * x["sem"]
    ci_r = lambda x: x["mean"] + 1.96 * x["sem"]

    # Get statistics
    stats = (
        df.groupby(groupby)["cadd_phred"]
        .agg(n="count", mean="mean", std=np.std, sem="sem")
        .assign(ci_l=ci_l, ci_r=ci_r)
    )

    return stats


def concat_transcript_data(df1, df2):
    return pd.concat([df1, df2], axis=0).sort_index()


def sort_region(df, **kwargs):
    kwargs.setdefault("categories", C.MAPS_CONSEQUENCES)
    kwargs.setdefault("labels", C.MAPS_LABELS)
    return (
        df.reset_index()
        .pipe(statistics_for_plots.sort_region_column, **kwargs)
        .set_index(["constraint", "region"])
        .sort_index()
    )


def main():
    # Read data
    df = read_cadd_annotations(C.CADD_ANNOTATED)

    # Synonymous, missense, and nonsense by region
    syn = df[df["csq"] == "synonymous_variant"].copy().pipe(get_cadd_stats)
    mis = df[df["csq"] == "missense_variant"].copy().pipe(get_cadd_stats)
    stop = df[df["csq"] == "stop_gained"].copy().pipe(get_cadd_stats)

    # Variants across the whole transcript
    whole_transcript = df.pipe(
        get_cadd_stats, groupby=["constraint", "csq"]
    ).rename_axis(["constraint", "region"])

    # Concatenate the regional and transcript-level data
    cat_sort = lambda x: concat_transcript_data(x, whole_transcript).pipe(sort_region)
    syn, mis, stop = [cat_sort(x) for x in [syn, mis, stop]]

    # Write to output, log T test statistics
    for df, csq in zip([syn, mis, stop], ["synonymous", "missense", "nonsense"]):
        df.to_csv(f"data/statistics/cadd_{csq}.tsv", sep="\t")

        logger.info(
            f"T test statistics for {csq} variants:\n"
            f"{statistics_for_plots.test_constrained_vs_unconstrained(df)}"
        )

    pass


if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)
    main()
