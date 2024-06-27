"""MAPS calculation."""

# Imports
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial as P
import pandas as pd

from src import setup_logger
from src import constants as C

# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def maps_model(mu, ps, n_obs, deg=1):
    """Fit a polynomial regression to PS vs mu for synonymous contexts.

    See notebooks/02_maps_model_choices.ipynb for further information.

    Args:
        deg (optional, int): Degree of the polynomial, passed to np.Polynomial.fit().
            Defaults to 1.
    """

    p = P.fit(mu, ps, w=n_obs, deg=deg)

    logger.info(f"Polynomial expression: {p}")

    return p


def maps(p, df, transform=False):
    """Calculate MAPS and confidence intervals.

    Args:
        polynomial (np.Polynomial): Polynomial object to predict PS.
        df (DataFrame): Dataframe containing the following columns:
            "mu" (Float): Mutability
            "ps" (Float): Proportion of singeltons
            "n_obs" (Int): Number of observed variants
        transform (optional, str): Type of transformation. Options:
            "log_ps": Log transformation of "ps"
            "sqrt_mu": Square root transformation of "mu"

    """
    if not transform in [False, "log_ps", "sqrt_mu"]:
        raise ValueError(
            "Incorrect value for transform. Expected one of: 'log_ps', 'sqrt_mu'."
        )

    if transform == "log_ps":
        df["ps_pred"] = np.exp(p(df["mu"]))

    if transform == "sqrt_mu":
        df["ps_pred"] = p(np.sqrt(df["mu"]))

    if not transform:
        df["ps_pred"] = p(df["mu"])

    df["maps"] = df["ps"] - df["ps_pred"]
    df["se"] = np.sqrt((df["ps"] * (1 - df["ps"])) / df["n_obs"])
    df["ci95"] = 1.96 * df["se"]

    return df


def weighted_average(x1, x2, w1, w2):
    return np.average([x1, x2], weights=[w1, w2], axis=0)


def combine_non_cpg_and_cpg_data(non_cpg, cpg):
    joint = pd.concat([non_cpg, cpg]).reset_index()

    weighted_mean = lambda x: np.average(x, weights=joint.loc[x.index, "n_obs"])

    joint = joint.groupby("csq").agg(
        mu = ("mu", weighted_mean),
        n_singletons = ("n_singletons", "sum"),
        n_obs = ("n_obs", "sum"),
        ps = ("ps", weighted_mean),
        ps_pred = ("ps_pred", weighted_mean),
        maps = ("maps", weighted_mean),
    )

    joint["se"] = np.sqrt((joint["ps"] * (1 - joint["ps"])) / joint["n_obs"])
    joint["ci95"] = 1.96 * joint["se"]

    logger.info(f"MAPS scores:\n{joint.maps}")
    
    return joint


def main():
    """Run as script."""

    # Read synonymous contexts
    syn = pd.read_csv(C.PS_SYN_CONTEXT, sep="\t").query("n_singletons > 0")

    # Split synonymous contexts into non-CpG and CpG
    syn_non = syn[syn.variant_type == "non-CpG"].copy()
    syn_cpg = syn[syn.variant_type == "CpG"].copy()

    # Get a polynomial object for predicting proportion of singletons
    non_p = maps_model(np.sqrt(syn_non.mu), syn_non.ps, syn_non.n_obs)
    cpg_p = maps_model(np.sqrt(syn_cpg.mu), syn_cpg.ps, syn_cpg.n_obs)

    # Get proportion of singletons for consequences in CpG and non-CpG contexts
    non = pd.read_csv(C.PS_REGIONS_NON_CPG, sep="\t")
    cpg = pd.read_csv(C.PS_REGIONS_CPG, sep="\t")

    # # Calculate MAPS
    non = maps(non_p, non, transform="sqrt_mu")
    cpg = maps(cpg_p, cpg, transform="sqrt_mu")

    # Combine non-CpG and CpG data
    joint = combine_non_cpg_and_cpg_data(non, cpg)

    # Write to output
    joint.to_csv(C.MAPS, sep="\t")

    return joint  #! Testing


if __name__ == "__main__":
    joint = main()
