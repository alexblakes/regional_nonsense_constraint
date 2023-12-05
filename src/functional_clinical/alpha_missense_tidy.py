""" 
Module docstring. 
"""

# Imports
from pathlib import Path

import pandas as pd

from src import constants as C
from src import setup_logger

# Module constants
_NAMES = [
    "chr",
    "pos",
    "ref",
    "alt",
    "alignment",
    "uniprot",
    "enst",
    "aa",
    "alpha_missense",
    "classification",
]
_USECOLS = ["chr", "pos", "ref", "alt", "enst", "alpha_missense"]

# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_alpha_missense(path):
    """Read Alpha Missense data."""

    logger.info("Reading Alpha Missense data.")

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        # nrows=100000, #! Testing
        names=_NAMES,
        usecols=_USECOLS,
    )

    logger.info(f"Alpha Missense annotations: {len(df)}")

    return df


def get_lowest_score_per_site(df):
    """Get the lowest AlphaMissense score per site and transcript ID."""

    logger.info("Getting lowest score per site.")

    df = (
        df.sort_values("alpha_missense", ascending=True)
        .drop_duplicates(["chr", "pos", "enst"], keep="first")
        .sort_values(["chr", "pos", "enst"])
        .rename(columns={"alpha_missense":"alpha_missense_min"})
    )

    logger.info(f"Remaining sites: {len(df)}")

    return df


def drop_version_numbers(transcript_id):
    return transcript_id.split(".")[0]


def tidy_alpha_missense_scores(df):
    """Drop irrelevant columns and remove version numbers from transcript ids."""

    df = df.drop(["ref", "alt"], axis=1)

    df["enst"] = df["enst"].apply(drop_version_numbers)

    logger.info(f"Duplicated by site / enst: {df.duplicated(['chr','pos','enst']).sum()}")
    logger.info(f"NaN values:\n{df.isna().sum()}")

    return df


def main():
    """Run the script."""

    df = (
        read_alpha_missense(C.ALPHA_MISSENSE)
        .pipe(get_lowest_score_per_site)
        .pipe(tidy_alpha_missense_scores)
    )

    # Write output
    logger.info("Writing to output.")
    df.to_csv(C.ALPHA_MISSENSE_TIDY, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    main()
