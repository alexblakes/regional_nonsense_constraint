"""Module docstring here."""
# Imports
from collections import defaultdict
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C

# Module constants
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_mutability_data(path):
    """Read gnomAD mutation rate data."""

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
    """Docstring."""

    return df.replace(_COMPLEMENT)


def reverse_complement_contexts(df, column="tri"):
    """Docstring."""
    df[column] = pd.Series(
        ["".join([_COMPLEMENT[N] for N in tri])[::-1] for tri in df[column]]
    )

    return df


def annotate_variant_types(df):
    """Docstring."""

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
    m7 = 

    return df


def main():
    mu = get_mutability_data(C.GNOMAD_NC_MUTABILITY)
    mu_rev = reverse_complement_alleles(mu).pipe(reverse_complement_contexts)
    mu = pd.concat([mu, mu_rev])
    # mu.to_csv()
    return mu  # TODO Testing


if __name__ == "__main__":
    main()

# Add variant type annotations
