"""
Find the expected number of variants for a given transcript / NMD region.
"""

# Imports
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C

# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_all_transcripts_for_cpg_variants(df, transcript_list):
    """Merge CpG variants with all transcripts.

    Fewer transcripts are represented among CpG variants, compared to non-CpG variants.
    """

    df = df.merge(transcript_list, how="right")

    # Fill any NaN values in csq and region with dummy values
    df["csq"] = df.csq.fillna("missense_variant")
    df["region"] = df.region.fillna("transcript")

    # Fill NaNs in variant_type
    df["variant_type"] = df["variant_type"].fillna("CpG")

    return df


def reindex_data(df, cpg=False):
    """Ensure all combinations of enst / region / csq are present.

    Args:
        cpg (bool, default=False): Processes CpGs if True, non-CpGs if False.
    """

    logger.info(f"Re-indexing variants: {df.variant_type.unique()}")
    logger.info(f"Number of unique transcripts: {df.enst.nunique()}")

    df = (
        df.set_index(["enst", "region", "csq"])
        .unstack("enst")
        .stack(dropna=False)
        .reset_index()
    )

    # Check that all combinations of enst, region, and csq are present
    logger.info(f"Number of rows: {len(df)}")
    logger.info(
        f"Number of combinations of enst/region/csq: {df.enst.nunique() * df.region.nunique() * df.csq.nunique()}"
    )

    # Fill NaN values for variant type
    if cpg:
        df["variant_type"] = df["variant_type"].fillna("CpG")

    elif not cpg:
        df["variant_type"] = df["variant_type"].fillna("non-CpG")

    # Logging
    logger.info(f"Region value counts:\n{df.region.value_counts()}")
    logger.info(f"CSQ value counts:\n{df.csq.value_counts()}")
    logger.info(f"Transcript counts: {df.enst.nunique()}")
    logger.info(f"Observed variants: {df.obs.count()}")

    return df


def main():
    """Run the script."""

    # ? Does the script account for duplicated PAR transcripts on X and Y

    # Read regional variant data
    regions = pd.read_csv(C.OBS_COUNTS_REGIONS_20, sep="\t")

    # Get list of transcripts
    all_transcripts = regions[["chr", "enst"]].drop_duplicates()
    logger.info(f"Total number of transcripts: {len(all_transcripts)}")

    # Split region variants into CpGs and non-CpGs
    non_cpg = regions[regions.variant_type == "non-CpG"].pipe(reindex_data, cpg=False)
    cpg = (
        regions[regions.variant_type == "CpG"]
        .pipe(get_all_transcripts_for_cpg_variants, all_transcripts)
        .pipe(reindex_data, cpg=True)
    )

    # Combine regions
    regions = pd.concat([non_cpg, cpg])

    # Write to output
    logger.info("Writing to output.")
    regions.to_csv(C.OBS_COUNTS_REGIONS_20_CLEAN, sep="\t", index=False)

    # return regions  #! Testing


if __name__ == "__main__":
    main()
