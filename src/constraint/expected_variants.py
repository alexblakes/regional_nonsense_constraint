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
def process_synonymous_variants(syn, cpg=False):
    """Retrieve synonymous variants in CpG or non-CpG contexts

    Args:
        cpg (bool, default=False): Processes CpGs if True, non-CpGs if False.
    """

    if not cpg:
        syn = syn[syn["variant_type"] != "CpG"]
        logger.info(f"Non-CpG synonymous contexts: {len(syn)}")

    if cpg:
        syn = syn[syn["variant_type"] == "CpG"].copy()

        # Saturated CpG contexts are dropped
        mask = syn["prop_obs"] == 1
        syn = syn[~mask]
        logger.info(f"Saturated CpG synonymous contexts: {mask.sum()}")
        logger.info(f"Remaining CpG synonymous contexts: {len(syn)}")

    return syn


def variant_expectation_model(df):
    """Variant expectation model.
    
    Predicts the proportion of possible variants observed, for a given mutability.
    """

    y = np.log(1 - df["prop_obs"])
    x = df["mu"]
    X = sm.tools.add_constant(x)
    w = df["pos"]

    model = sm.WLS(y, X, weights=w)
    results = model.fit()
    logger.info(f"\n{results.summary(slim=True)}")

    return results

def reindex_non_cpgs(df):
    """Ensure all combinations of enst / region / csq are present."""
    
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
    logger.info(f"Number of combinations of enst/region/csq: {df.enst.nunique() * df.region.nunique() * df.csq.nunique()}")

    # Fill NaN values for variant type
    df["variant_type"] = df["variant_type"].fillna("non-CpG")

    return df


def main():
    """Run the script."""

    # Read synonymous counts
    syn = pd.read_csv(C.OBS_COUNTS_SYN_20, sep="\t")
    
    # Get proportion observed
    syn["prop_obs"] = syn["obs"] / syn["pos"]

    logger.info(f"Total synonymous contexts: {len(syn)}")

    # Split CpG and non-CpG variant contexts
    syn_non = process_synonymous_variants(syn, cpg=False)
    syn_cpg = process_synonymous_variants(syn, cpg=True)

    non_cpg_results = variant_expectation_model(syn_non)
    cpg_results = variant_expectation_model(syn_cpg)

    # syn_non = reindex_non_cpgs(syn_non)

    return syn_non  #! Testing


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
