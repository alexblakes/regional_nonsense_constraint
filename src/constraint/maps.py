"""MAPS calculation."""

# Imports
from pathlib import Path

import numpy as np
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
        mu (float): Mutability.
        ps (float): Proportion of singletons.
        n_obs (int): Number of variants observed.
        deg (optional, int): Degree of the polynomial, passed to np.polyfit().
    """

    z = np.polyfit(mu, ps, w=n_obs, deg=deg)
    p = np.poly1d(z)

    return p


def maps(poly1d, df, transform=False):
    """Calculate MAPS and confidence intervals.

    Args:
        poly1d (np.poly1d): Polynomial class to predict PS.
        df (DataFrame): Dataframe containing the following columns:
            "mu" (Float): Mutability
            "ps" (Float): Proportion of singeltons
            "n_obs" (Int): Number of observed variants
        transform (optional, str): Type of transformation. Options:
            "log_ps": Log transformation of "ps"
            "sqrt_mu": Square root transformation of "mu"
            
    """

    if transform == "log_ps":
        df["ps_pred"] = np.exp(poly1d(df["mu"]))

    if transform == "sqrt_mu":
        df["ps_pred"] = poly1d(np.sqrt(df["mu"]))

    if not transform:
        df["ps_pred"] = poly1d(df["mu"])
    
    df["maps"] = df["ps"] - df["ps_pred"]
    df["se"] = np.sqrt((df["ps"] * (1 - df["ps"])) / df["n_obs"])
    df["ci95"] = 1.96 * df["se"]

    return df


def main():
    """Run as script."""

    syn = (
        pd.read_csv(C.PS_SYN_CONTEXT, sep="\t").query("n_singletons > 0")
        # .query("n_obs >= 1000")
        # .query("mu < (9 * 10 ** -8)")
    )

    # Split into CpG and non-CpG contexts
    cpg_p = syn[syn.variant_type == "CpG"].pipe(maps_model, deg=2)
    non_p = syn[syn.variant_type == "non-CpG"].pipe(maps_model, deg=1)

    # Calculate MAPS separately for CpG and non-CpG variants
    cpg = pd.read_csv(C.PS_REGIONS_CPG, sep="\t").pipe(maps, cpg_p)
    non = pd.read_csv(C.PS_REGIONS_NON_CPG, sep="\t").pipe(maps, non_p)

    p = maps_model(syn, deg=2)
    joint = pd.read_csv(C.PS_REGIONS, sep="\t").pipe(maps, p)

    return cpg, non, joint  #! Testing

if __name__ == "__main__":
    cpg, non, joint = main()