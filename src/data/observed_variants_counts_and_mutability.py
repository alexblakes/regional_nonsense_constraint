"""
Get variant counts and mean mutability for synonymous contexts and NMD regions.

These statistics are created for rare synonymous variants, and variants in each NMD 
region per transcript.
"""

# Imports
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C


# Module constants
# "categorical" dtype is not parsed from a defaultdict.
_DATATYPES = {
    "chr": "category",
    "pos": "int32",
    "ref": "category",
    "alt": "category",
    "csq": "category",
    "enst": "category",
    "tri": "category",
    "nmd": "category",
    "ac": "Int32",
    "an": "Int32",
    "obs": "bool",
    "variant_type": "category",
    "lvl": "category",
    "mu": "float32",
    "median_coverage": "Int16",
}

# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_variant_annotations(path):
    """Read variant annotation data."""

    logger.info("Reading variant annotation data.")

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=_DATATYPES,
        # nrows=10000,  #! Testing
    )

    return df


def get_rare_synonymous_variants(df):
    """Get rare synonymous variants for variant expectation model."""

    logger.info("Filtering for rare synonymous variants.")

    syn = df[df["csq"] == "synonymous_variant"]

    # Mask contexts in which a synonymous variant is generally not possible.
    # Synonymous variants in these contexts can only occur at exon-intron junctions.
    m1 = (syn.tri == "AGT") & ((syn.alt == "C") | (syn.alt == "T"))
    m2 = (syn.tri == "AAT") & ((syn.alt == "C") | (syn.alt == "T"))
    m3 = (syn.tri == "ACT") & ((syn.alt == "G") | (syn.alt == "A"))
    m4 = (syn.tri == "ATT") & ((syn.alt == "G") | (syn.alt == "A"))

    # Mask rare variants
    m5 = syn["ac"].isna()  # i.e. ac == 0
    m6 = (syn["ac"] / syn["an"]) < 0.001  # i.e. AF < 0.1%

    # Get qualifying synonymous variants
    syn = syn[~(m1 | m2 | m3 | m4) & (m5 | m6)]

    # Logging
    logger.info(f"Synonymous variants at splice junctions: {(m1 | m2 | m3 | m4).sum()}")
    logger.info(f"Synonymous variants not observed: {m5.sum()}")
    logger.info(f"Observed rare (AF < 0.1%) synonymous variants: {m6.sum()}")
    logger.info(f"Qualifying synonymous variants: {len(syn)}")

    return syn


def agg_counts_and_mutability(df, grouping_columns):
    """Get variant counts and mean mutability for a given grouping."""

    df = (
        df.groupby(grouping_columns, observed=True)
        .agg({"obs": "sum", "pos": "count", "mu": "mean"})
        .reset_index()
    )

    return df


def main():
    """Run the script."""

    df = get_variant_annotations(C.ALL_VARIANTS_MERGED_ANNOTATIONS)
    # ? We could filter for coverage at this point

    # Synonymous variants
    syn = get_rare_synonymous_variants(df).pipe(
        agg_counts_and_mutability, ["tri", "ref", "alt", "variant_type", "lvl"]
    )

    # Regions
    # We need to group by chr for PAR transcripts on chrX & chrY
    transcript = agg_counts_and_mutability(
        df, ["chr", "enst", "csq", "variant_type"]
    ).assign(nmd="transcript")
    nmd = agg_counts_and_mutability(df, ["chr", "enst", "csq", "variant_type", "nmd"])
    regions = pd.concat([nmd, transcript])

    # Logging
    logger.info(f"Synonymous contexts annotated: {len(syn)}")
    logger.info(f"enst/csq/cpg annotations: {len(transcript)}")
    logger.info(f"nmd/enst/csq/cpg annotations: {len(nmd)}")

    # Write to output
    syn.to_csv(C.OBSERVED_VARIANTS_COUNTS_SYN, sep="\t", index=False)
    regions.to_csv(C.OBSERVED_VARIANTS_COUNTS_REGION, sep="\t", index=False)

    # return regions  #! Testing


if __name__ == "__main__":
    main()
