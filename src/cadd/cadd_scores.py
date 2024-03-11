"""Annotate sites with CADD scores."""

# Imports
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C

# Module constants
_CADD_HEADER = "chr pos ref alt cadd_raw cad_phred".split()
_CADD_USECOLS = "chr pos ref alt cad_phred".split()

# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_cadd(path):
    """Read CADD data."""

    cadd = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=_CADD_HEADER,
        usecols=_CADD_USECOLS,
        nrows=10,
    )

    return cadd

def main():
    """Run as script."""

    cadd = read_cadd("data/raw/whole_genome_SNVs.tsv.gz")

    return cadd

if __name__ == "__main__":
    main()
