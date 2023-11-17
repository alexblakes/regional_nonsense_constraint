"""
Find the expected number of variants for a given transcript / NMD region.
"""

# Imports
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

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

    regions = pd.read_csv(C.OBS_COUNTS_REGIONS_20, sep="\t")

    # Get list of transcripts
    all_transcripts = pd.Series(regions["enst"].unique(), name="enst")
    logger.info(f"Total number of transcripts: {len(all_transcripts)}")

    # Split region variants into CpGs and non-CpGs
    non_cpg = regions[regions.variant_type == "non-CpG"].pipe(reindex_data, cpg=False)
    cpg = (
        regions[regions.variant_type == "CpG"]
        .pipe(get_all_transcripts_for_cpg_variants, all_transcripts)
        .pipe(reindex_data, cpg=True)
    )

    return cpg  #! Testing


if __name__ == "__main__":
    main()


# # ## Calculate the expected proportion of variants per transcript and region


# # ### Non-CpGs


# # Predict the proportion observed with the non-CpG model
# non_cpg["prop_exp"] = 1 - np.exp(non_cpg_results.predict(sm.tools.add_constant(non_cpg["mu"])))

# # Calculate the number of variants expected
# non_cpg["n_exp"] = (non_cpg["prop_exp"] * non_cpg["n_pos"]).pipe(np.round, 3)


# # ### CpGs


# # Predict the proportion observed with the CpG model
# cpg["prop_exp"] = 1 - np.exp(cpg_results.predict(np.exp(sm.tools.add_constant(cpg["mu"]))))

# # Calculate the number of variants expected
# cpg["n_exp"] = (cpg["prop_exp"] * cpg["n_pos"]).pipe(np.round, 3)


# # ### Combine CpG and non-CpG variants


# # Give CpG and non-CpG variants a consistent index
# non_cpg = non_cpg.set_index(["enst","region","csq"])
# cpg = cpg.set_index(["enst","region","csq"])

# # Combine non-CpG and CpG data into the essential summary statistics
# n_pos = non_cpg["n_pos"].fillna(0) + cpg["n_pos"].fillna(0)
# n_obs = non_cpg["n_obs"].fillna(0) + cpg["n_obs"].fillna(0)
# n_exp = non_cpg["n_exp"].fillna(0) + cpg["n_exp"].fillna(0)
# oe = (n_obs / n_exp).rename("oe")
# prop_obs = (n_obs / n_pos).rename("prop_obs")
# prop_exp = (n_exp / n_pos).rename("prop_exp")

# # Calculate the total mutability for each region
# # In each dataframe, "mu" is the mean mutability for contexts in a region
# mu_non_cpg = (non_cpg["mu"] * non_cpg["n_pos"]).fillna(0)
# mu_cpg = (cpg["mu"] * cpg["n_pos"]).fillna(0)
# mu = (mu_non_cpg + mu_cpg).rename("mu")

# # Create a summary dataframe
# df = pd.concat([mu, n_pos, n_obs, n_exp, oe, prop_obs, prop_exp], axis=1).reset_index(drop=False)

# # How many regions are missing variants?
# oe.isna().sum()


# #
# # Note the many sites with NaN values. These are regions where 0 variants are possible. We have kept them for now; they could be dropped later.


# # ## Write to output


# df.to_csv("../outputs/expected_variants_all_regions.tsv", sep="\t", index=False)
