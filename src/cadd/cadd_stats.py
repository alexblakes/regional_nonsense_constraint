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
_CSQS = "synonymous missense nonsense".split()
_REGION_LABELS = [
    "Synonymous (Whole CDS)",
    "Missense (Whole CDS)",
    "Nonsense (Whole CDS)",
    "NMD target",
    "Start proximal",
    "Long exon",
    "Distal",
][
    ::-1
]  # Reversed for plotting

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
            # nrows=10000000,
        )
        .dropna(subset="constraint")
        .replace({"region": "distal_nmd"}, value="distal")
    )


def get_cadd_stats(df, groupby=["constraint", "region"]):
    """Get summary stats of CADD scores by region and constraint."""

    # Lambda functions for 95% confidence intervals
    ci_95 = lambda x: 1.96 * x["sem"]

    # Get statistics
    stats = (
        df.groupby(groupby)["cadd_phred"]
        .agg(n="count", mean="mean", std=np.std, sem="sem")
        .assign(ci_95=ci_95)
    )

    return stats


def concat_transcript_data(df1, df2):
    return pd.concat([df1, df2], axis=0).sort_index()


def tidy_region_names(df, **kwargs):
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
    cat_sort = lambda x: concat_transcript_data(x, whole_transcript).pipe(
        tidy_region_names, labels=_REGION_LABELS
    )
    syn, mis, stop = [cat_sort(x) for x in [syn, mis, stop]]
    logger.debug(f"{syn}")
    logger.debug(f"{mis}")
    logger.debug(f"{stop}")

    # Outputs
    for df, path, csq in zip(
        [syn, mis, stop], [C.STATS_CADD_SYN, C.STATS_CADD_MIS, C.STATS_CADD_NON], _CSQS
    ):
        # Write to output
        df.to_csv(path, sep="\t")

        # Log T test statistics
        logger.info(
            f"T test statistics for {csq} variants:\n"
            f"{statistics_for_plots.test_constrained_vs_unconstrained(df)}"
        )

    pass


if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)
    main()
