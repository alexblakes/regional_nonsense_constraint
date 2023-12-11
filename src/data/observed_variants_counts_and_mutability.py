"""
Get variant counts and mean mutability for synonymous contexts and NMD regions.

These statistics are created for rare synonymous variants, and variants in each NMD 
region per transcript.
"""

# Imports
import argparse
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


def filter_covered_sites(df, coverage):
    """Filter for sites with minimum coverage requirement."""

    if not type(coverage) == int:
        raise ValueError("Coverage must be integer.")

    logger.info(f"Filtering for sites with coverage >= {coverage}.")

    def covered_sites_logging(df):
        logger.info(f"Qualifying variants: {len(df)}")
        logger.info(f"Observed variants: {df.obs.sum()}")
        logger.info(f"Qualifying variants by consequence:\n{df.csq.value_counts()}")
        logger.info(
            f"Qualifying variants by observed / consequence: {df.groupby('obs').csq.value_counts()}"
        )

    # Some sites have NaN values for coverage.
    # If coverage == 0, do not exclude even these.
    if coverage == 0:
        logger.info(
            f"Variants where coverage is NaN: {df.median_coverage.isna().sum()}"
        )
        covered_sites_logging(df)
        return df

    df = df[df.median_coverage >= coverage]
    covered_sites_logging(df)

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


def parse_args():
    """Parse command line arguments."""

    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument(
        "-c",
        "--coverage",
        type=int,
        default=[0],
        nargs="*",
        help="Minimum coverage of sites to include. Accepts multiple integer values.",
    )

    return parser.parse_args()


def main():
    """Run the script."""

    df = get_variant_annotations(C.ALL_VARIANTS_MERGED_ANNOTATIONS).rename(
        columns={"nmd": "region"}
    )

    # Filter for coverage
    coverage = parse_args().coverage

    for cov in coverage:
        # Filter variants by coverage
        df_min_coverage = filter_covered_sites(df, cov)

        # Get rare synonymous variants for expectation model
        syn = get_rare_synonymous_variants(df_min_coverage).pipe(
            agg_counts_and_mutability, ["tri", "ref", "alt", "variant_type", "lvl"]
        )

        # Get variants by region for constraint calculations.
        # ? Should we limit to rare variants only?
        # * We need to group by chr for PAR transcripts on chrX & chrY
        transcript = agg_counts_and_mutability(
            df_min_coverage, ["chr", "enst", "csq", "variant_type"]
        ).assign(region="transcript")

        nmd = agg_counts_and_mutability(
            df_min_coverage, ["chr", "enst", "csq", "variant_type", "region"]
        )
        regions = pd.concat([nmd, transcript])

        # Logging
        logger.info(f"Synonymous contexts annotated: {len(syn)}")
        logger.info(f"Counts grouped by enst/csq/cpg: {len(transcript)}")
        logger.info(f"Counts grouped by nmd/enst/csq/cpg: {len(nmd)}")

        # Write to output
        logger.info("Writing to output.\n")
        syn.to_csv(f"{C._OBS_COUNTS_SYN}{str(cov)}.tsv", sep="\t", index=False)
        regions.to_csv(f"{C._OBS_COUNTS_REGIONS}{str(cov)}.tsv", sep="\t", index=False)

    # return regions  #! Testing


if __name__ == "__main__":
    main()
