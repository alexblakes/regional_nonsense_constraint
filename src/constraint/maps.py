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
def fit_maps_model(syn):
    """Fit a second degree polynomial to PS vs sqrt(mu).

    See notebooks/02_maps_model_choices.ipynb for further information.
    """

    x = np.sqrt(syn.mu)
    y = syn.ps
    w = syn.n_obs

    z = np.polyfit(x, y, deg=2, w=w)
    p = np.poly1d(z)

    return p


def maps(df, poly1d):
    """Calculate MAPS and confidence intervals.

    Args:
        poly1d (np.poly1d): Polynomial class to predict PS.
    """

    df["ps_pred"] = poly1d(np.sqrt(df.mu))
    df["maps"] = (df["ps"] - df["ps_pred"]).pipe(np.round, 3)
    df["se"] = np.sqrt((df["ps"] * (1 - df["ps"])) / df["n_obs"])
    df["ci95"] = 1.96 * df["se"]

    return df


def main():
    """Run as script."""

    df = pd.read_csv(C.PS_REGIONS, sep="\t")
    p = pd.read_csv(C.PS_SYN_CONTEXT, sep="\t").pipe(fit_maps_model)

    df = maps(df, p)

    return df  #! Testing
