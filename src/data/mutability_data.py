"""Tidy the gnomAD non-coding constraint mutability data."""

# Imports
from pathlib import Path

import numpy as np
import pandas as pd

from src import setup_logger
from src import constants as C

# Module constants
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_mutability_data(path):
    """Read mutability data."""

    logger.info(f"Reading mutability data from {C.GNOMAD_NC_MUTABILITY}")

    mu = pd.read_csv(
        path,
        sep="\t",
        names=[
            "tri",
            "ref",
            "alt",
            "lvl",
            "pos",
            "obs",
            "po",
            "mu",
            "ppo",
        ],
        header=0,
        usecols=["tri", "ref", "alt", "lvl", "mu"],
    )

    return mu


def reverse_complement_alleles(df):
    """Replace alleles in the REF and ALT columns with their reverse-complement."""

    logger.info("Getting reverse complement of REF / ALT alleles")    
    return df.replace(_COMPLEMENT)


def reverse_complement_contexts(df, column="tri"):
    """Reverse complement trinucleotide contexts."""
    
    logger.info("Getting reverse complement of trinucleotide contexts.")    
    
    df[column] = pd.Series(
        ["".join([_COMPLEMENT[N] for N in tri])[::-1] for tri in df[column]]
    )

    return df


def annotate_variant_types(df):
    """Annotate CpG and non-CpG variants."""

    logger.info("Annotating CpG variants.")

    # Masks for CpGs transitions
    ## "Forward" strand
    m1 = df.tri.str[1] == "C"
    m2 = df.tri.str[2] == "G"
    m3 = df.alt == "T"

    ## "Reverse" strand
    m4 = df.tri.str[0] == "C"
    m5 = df.tri.str[1] == "G"
    m6 = df.alt == "A"

    ## Combined
    m7 = (m1 & m2 & m3) | (m4 & m5 & m6)

    df["variant_type"] = np.where(m7, "CpG", "non-CpG")

    return df


def main():
    """Run the script."""

    mu = get_mutability_data(C.GNOMAD_NC_MUTABILITY)
    mu_rev = reverse_complement_alleles(mu).pipe(reverse_complement_contexts)
    mu = pd.concat([mu, mu_rev]).pipe(annotate_variant_types)

    # Logging
    logger.info(f"Number of variant contexts: {len(mu)}")
    logger.info(f"Check for duplicates: {mu.duplicated().sum()}")
    logger.info(f"Number of CpG contexts: {(mu.variant_type == 'CpG').sum()}")
    logger.info("Writing to TSV.")

    mu.to_csv(C.GNOMAD_NC_MUTABILITY_TIDY, sep="\t", index=False)

    return mu  # TODO Testing


if __name__ == "__main__":
    main()
