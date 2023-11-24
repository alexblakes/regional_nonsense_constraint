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
def maps_model(syn, deg):
    """Fit a second degree polynomial to PS vs mu for synonymous contexts.

    See notebooks/02_maps_model_choices.ipynb for further information.
    """

    x = syn.mu
    y = syn.ps
    w = syn.n_obs

    z = np.polyfit(x, y, w=w, deg=deg)
    p = np.poly1d(z)

    return p


# def fit_maps_model(syn):
#     """Fit a second degree polynomial to PS vs sqrt(mu).

#     See notebooks/02_maps_model_choices.ipynb for further information.
#     """

#     x = np.sqrt(syn.mu)
#     y = syn.ps
#     w = syn.n_obs

#     z = np.polyfit(x, y, deg=2, w=w)
#     p = np.poly1d(z)

#     return p


def maps(df, poly1d):
    """Calculate MAPS and confidence intervals.

    Args:
        poly1d (np.poly1d): Polynomial class to predict PS.
    """

    df["ps_pred"] = poly1d(df.mu)
    df["maps"] = (df["ps"] - df["ps_pred"]).pipe(np.round, 3)
    df["se"] = np.sqrt((df["ps"] * (1 - df["ps"])) / df["n_obs"])
    df["ci95"] = 1.96 * df["se"]

    return df


def main():
    """Run as script."""

    syn = (
        pd.read_csv(C.PS_SYN_CONTEXT, sep="\t")
        .query("n_singletons > 0")
        .query("n_obs >= 1000")
    )

    # Split into CpG and non-CpG contexts
    cpg_p = syn[syn.variant_type == "CpG"].pipe(maps_model, deg=2)
    non_p = syn[syn.variant_type == "non-CpG"].pipe(maps_model, deg=1)

    # Calculate MAPS separately for CpG and non-CpG variants
    cpg = pd.read_csv(C.PS_REGIONS_CPG, sep="\t").pipe(maps, cpg_p)
    non = pd.read_csv(C.PS_REGIONS_NON_CPG, sep="\t").pipe(maps, non_p)

    return cpg, non  #! Testing
