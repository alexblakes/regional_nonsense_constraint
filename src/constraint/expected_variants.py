""" Find the expected number of variants in all transcripts / NMD regions."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C

_FILE_IN = "data/interim/observed_variants_counts_regions_cov_20.tsv"
_FILE_OUT = "data/final/expected_variants_all_regions.tsv"
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_NAMES = "enst csq region n_obs n_pos mu".split()
_SCALE = 1.015 * (10**-7)  # Scaling factor for per-generation mutation rate

logger = logging.getLogger(__name__)


def read_variant_counts(path):
    return pd.read_csv(path, sep="\t", header=0, names=_NAMES)


def assign_statistics(df):
    return (
        df.assign(prop_exp=lambda x: x["mu"])
        .assign(mu=lambda x: x["mu"] * _SCALE)
        .assign(n_exp=lambda x: x["n_pos"] * x["prop_exp"])
        .assign(oe=lambda x: x["n_obs"] / x["n_exp"])
        .assign(prop_obs=lambda x: x["n_obs"] / x["n_pos"])
    )


def reorder_data(df):
    return df[
        [
            "enst",
            "region",
            "csq",
            "n_pos",
            "n_obs",
            "n_exp",
            "oe",
            "prop_obs",
            "prop_exp",
            "mu",
        ]
    ]


def main():
    """Run as script."""

    df = read_variant_counts(_FILE_IN).pipe(assign_statistics).pipe(reorder_data)

    # Logging
    logger.info(f"Number of entries: {len(df)}")
    logger.info(f"NaN values:\n{df.isna().sum()}")
    logger.info(f"Unique transcripts: {df.enst.nunique()}")
    logger.info(f"Region value counts:\n{df.groupby('region').csq.value_counts()}")

    # Write out
    df.to_csv(_FILE_OUT, sep="\t", index=False)

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
