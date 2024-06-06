""" Find the expected number of variants in all transcripts / NMD regions."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C

_FILE_IN = "data/interim/observed_variants_counts_regions_cov_20.tsv"
_FILE_OUT = "data/final/expected_variants_all_regions.tsv"
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

logger = logging.getLogger(__name__)

def read_variant_counts(path):
    return pd.read_csv(path, sep="\t")

def main():
    """Run as script."""

    df = read_variant_counts(_FILE_IN)

    return df #! Testing


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
