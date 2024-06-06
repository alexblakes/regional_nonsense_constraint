"""Find constrained and unconstrained regions."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C

_FILE_IN = "data/final/regional_constraint_stats.tsv"
_FILE_OUT = "data/final/regional_nonsense_constraint.tsv"
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_SYN_Z = -1
_OE = 0.35
_FDR_P = 0.05
_P = 0.05
_OBS = 1

logger = logging.getLogger(__name__)


def get_synonymous_z_scores(df):
    """Add synonymous z scores to nonsense constraint data.

    Return nonsense constraint data only.
    """

    # Get synonymous and nonsense variants
    syn = df[df.csq == "synonymous_variant"].copy().rename(columns={"z": "syn_z"})
    stop = df[df.csq == "stop_gained"].copy()

    # Annotate nonsense constraint data with synonymous z scores
    stop = stop.merge(syn[["enst", "region", "syn_z"]], how="inner")

    return stop


def assign_constraint_label(df):
    """Label constrained and unconstrained regions."""

    # Filtering masks
    m1 = df["oe"] < _OE
    m2 = df["syn_z"] >= _SYN_Z
    m3 = df["fdr_p"] < _FDR_P

    m4 = df["p"] >= _P
    m5 = df["n_obs"] >= _OBS

    # Assign constraint label
    df.loc[m1 & m2 & m3, "constraint"] = "constrained"
    df.loc[m4 & m5, "constraint"] = "unconstrained"

    # Logging
    logger.info(
        f"Constrained region value counts:\n{df.groupby('region').constraint.value_counts(dropna=False)}"
    )

    return df


def main():
    """Run as script."""

    # Read the data
    df = (
        pd.read_csv(_FILE_IN, sep="\t")
        .pipe(get_synonymous_z_scores)
        .pipe(assign_constraint_label)
    )

    # Write to output
    df.to_csv(_FILE_OUT, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
