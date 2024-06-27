"""Get the proportion of singletons for synonymous variant contexts."""

# Imports
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C

# Module constants
_COVERAGE = 20


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_variant_annotations(path):
    """Get merged annotations for all variants."""

    logger.info("Getting variant annotations.")

    df = (
        pd.read_csv(
            path,
            sep="\t",
            # nrows=1000000,
        )
        .rename(columns={"nmd": "region"})
        .query("obs == True")
        .query(f"median_coverage >= {_COVERAGE}")
    )

    # Sanity checks
    assert all(x > 0 for x in df.ac), "0s or NaNs present in df.ac"

    # Logging
    logger.info(f"Observed variants at {_COVERAGE}x coverage: {len(df)}")
    logger.info(f"Missing values:\n{df.isna().sum()}")

    return df


def get_valid_synonymous_variants(df):
    """Get valid synonymous variants for MAPS model."""

    logger.info("Filtering for valid synonymous variants.")

    syn = df[df["csq"] == "synonymous_variant"]

    logger.info(f"Synonymous variants before filtering: {len(syn)}")

    # Mask contexts in which a synonymous variant is generally not possible.
    # Synonymous variants in these contexts can only occur at exon-intron junctions.
    m1 = (syn.tri == "AGT") & ((syn.alt == "C") | (syn.alt == "T"))
    m2 = (syn.tri == "AAT") & ((syn.alt == "C") | (syn.alt == "T"))
    m3 = (syn.tri == "ACT") & ((syn.alt == "G") | (syn.alt == "A"))
    m4 = (syn.tri == "ATT") & ((syn.alt == "G") | (syn.alt == "A"))

    syn = syn[~(m1 | m2 | m3 | m4)]

    logger.info(f"Synonymous variants after filtering: {len(syn)}")

    return syn


def get_ps(dfg):
    """Get mean mutability and proportion of singletons."""

    logger.info("Getting proportion of singletons.")

    mu = dfg["mu"].mean()
    ns = dfg["ac"].apply(lambda x: (x == 1).sum()).rename("n_singletons")
    no = dfg["ac"].count().rename("n_obs")

    ps = pd.concat([mu, ns, no], axis=1).reset_index()

    ps["ps"] = ps["n_singletons"] / ps["n_obs"]

    return ps


def main():
    """Run as script."""

    df = (
        get_variant_annotations(C.ALL_VARIANTS_MERGED_ANNOTATIONS)
        .pipe(get_valid_synonymous_variants)
        .groupby(["tri", "ref", "alt", "variant_type", "lvl"])
        .pipe(get_ps)
    )

    # Write to output
    logger.info("Writing to output.")
    df.to_csv(C.PS_SYN_CONTEXT, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    main()
