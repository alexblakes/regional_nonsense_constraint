"""Module docstring here."""
# Imports
from collections import defaultdict
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C

# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_mutability_data(path):
    """Read gnomAD mutation rate data."""

    context	ref	alt	methylation_level	possible	observed	proportion_observed	mu	fitted_po


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



def main():
    df = get_mutability_data(C.GNOMAD_NC_MUTABILITY)

    return df # TODO Testing

if __name__ == "__main__":
    main()

# # Simplify variant type annotations
# mu = mu

# # Mutation rates are only available for 32 codons. We need to reverse-complement for the remainder.
# complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

# # Replace ref and alt alleles
# _mu = mu.copy().replace(complement)

# # Reverse-complement trinucleotide contexts
# _mu["tri"] = pd.Series(["".join([complement[y] for y in x])[::-1] for x in mu.tri])

# # Merge original and reverse-complemented data
# mu = pd.concat([mu, _mu])