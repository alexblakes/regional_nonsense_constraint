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
_REGION_CATEGORIES = "whole_cds nmd_target start_proximal long_exon distal".split()[
    ::-1
]  # Reversed for plotting
_REGION_LABELS = ["Whole CDS", "NMD target", "Start proximal", "Long exon", "Distal"][
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
            # nrows=2000000,
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


def concat_cds_data(region, cds, csq):
    # Subset CDS data by consequence
    cds = (
        cds.xs(csq, level="region", drop_level=False)
        .reset_index("region")
        .assign(region="whole_cds")
        .set_index("region", append=True)
    )

    return pd.concat([region, cds], axis=0).sort_index()


def tidy_region_names(df, **kwargs):
    kwargs.setdefault("categories", _REGION_CATEGORIES)
    kwargs.setdefault("labels", _REGION_LABELS)

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
    non = df[df["csq"] == "stop_gained"].copy().pipe(get_cadd_stats)

    # Variants across the whole transcript
    cds = df.pipe(get_cadd_stats, groupby=["constraint", "csq"]).rename_axis(
        ["constraint", "region"]
    )

    # Concatenate the regional and transcript-level data
    syn = concat_cds_data(syn, cds, "synonymous_variant").pipe(tidy_region_names)
    mis = concat_cds_data(mis, cds, "missense_variant").pipe(tidy_region_names)
    non = concat_cds_data(non, cds, "stop_gained").pipe(tidy_region_names)

    # Outputs
    for df, path, csq in zip(
        [syn, mis, non], [C.STATS_CADD_SYN, C.STATS_CADD_MIS, C.STATS_CADD_NON], _CSQS
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
    logger = src.setup_logger(_LOGFILE)
    main()
