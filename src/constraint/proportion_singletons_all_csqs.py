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
def main():
    """Run as script."""

    # Read variant annotations
    df = ps.get_variant_annotations(C.ALL_VARIANTS_MERGED_ANNOTATIONS)

    # Get PS per consequence
    ps_csq = df.groupby("csq").pipe(ps.get_ps)

    # For nonsense variants, get PS per region
    ps_region = (
        df[df["csq"] == "stop_gained"]
        .groupby("region")
        .pipe(ps.get_ps)
        .rename(columns={"region": "csq"})
    )

    # Combine datasets
    df = pd.concat([ps_csq, ps_region])

    # Write to output
    df.to_csv(C.PS_REGIONS, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    main()
