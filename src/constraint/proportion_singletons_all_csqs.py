"""Get the proportion of singletons for all variant classes."""

# Imports
from pathlib import Path

import pandas as pd

from src.constraint import proportion_singletons_syn_contexts as ps
from src import setup_logger
from src import constants as C

# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_ps_per_consequence(df):
    return df.groupby("csq").pipe(ps.get_ps)


def get_ps_per_region(df):
    return (
        df[df["csq"] == "stop_gained"]
        .groupby("region")
        .pipe(ps.get_ps)
        .rename(columns={"region": "csq"})
    )


def main():
    """Run as script."""

    # Read variant annotations
    df = ps.get_variant_annotations(C.ALL_VARIANTS_MERGED_ANNOTATIONS)

    # Split by variant type
    cpg = df[df.variant_type == "CpG"]
    non = df[df.variant_type == "non-CpG"]

    # Get PS per consequence
    cpg_csq = get_ps_per_consequence(cpg)
    non_csq = get_ps_per_consequence(non)

    # For nonsense variants, get PS per region
    cpg_region = get_ps_per_region(cpg)
    non_region = get_ps_per_region(non)

    # Combine datasets
    cpg = pd.concat([cpg_csq, cpg_region])
    non = pd.concat([non_csq, non_region])

    # Write to output
    logger.info("Writing to output.")
    cpg.to_csv(C.PS_REGIONS_CPG, sep="\t", index=False)
    non.to_csv(C.PS_REGIONS_NON_CPG, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    main()
