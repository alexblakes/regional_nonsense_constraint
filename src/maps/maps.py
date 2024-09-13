"""MAPS calculation."""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/interim/maps_snvs_filtered.tsv"
_FILE_OUT = "data/interim/maps.tsv"
_USECOLS = ["csq", "region", "ac", "mu_scaled"]
_DTYPES = {
    "csq": "category",
    "region": "category",
    "ac": np.float32,
    "mu_scaled": np.float32,
}

logger = logging.getLogger(__name__)


def read_data(path=_FILE_IN):
    df = pd.read_csv(
        path,
        sep="\t",
        na_values=".",
        header=0,
        usecols=_USECOLS,
        dtype=_DTYPES,
        # nrows=10000,
    )

    return df


def get_ps(df, grouping="mu_scaled", mu="mu_scaled"):
    """Get proportion of singletons."""

    grouped = df.groupby(grouping)

    ps = (
        grouped.agg(
            n_variants=("ac", "count"),
            n_singletons=("ac", lambda x: (x == 1).sum()),
            mu=(mu, "mean"),
        )
        .assign(ps=lambda x: x["n_singletons"] / x["n_variants"])
        .reset_index(drop=False)
    )

    return ps


def fit_lm(df, weighted=True, **kwargs):
    """Fit a polynomial to proportion singletons vs mutation rate."""

    kwargs.setdefault("deg", 1)

    x = df["mu"]
    y = df["ps"]

    if weighted:
        kwargs.setdefault("w", df["n_variants"])

    z = np.polyfit(x, y, **kwargs)

    return np.poly1d(z)


def pred_ps(df, poly1d):
    """Find the expected proportion of singletons for a given mutation rate."""

    if not poly1d:
        # For modelling synonymous variants
        poly1d = fit_lm(df)

    df["pred_ps"] = poly1d(df["mu"])

    return df


def get_maps(df):
    return df.assign(maps=lambda x: x["ps"] - x["pred_ps"])


def get_confidence_interval(df):
    df["ci95"] = 1.96 * np.sqrt((df["ps"] * (1 - df["ps"])) / df["n_variants"])
    return df


def main():
    """Run as script."""

    df = read_data()

    # Subset to synonymous variants only
    syn = df[df["csq"] == "synonymous_variant"].copy()

    logger.info(f"Scaled mutation rate levels: {syn['mu_scaled'].nunique()}")

    # Fit a linear model of PS ~ mu for synonymous variants
    # See related notebook for details
    p = syn.pipe(get_ps).pipe(fit_lm, weighted=True)

    # Calculate MAPS at the transcript level and regional level
    transcript = (
        get_ps(df, grouping="csq")
        .pipe(pred_ps, poly1d=p)
        .pipe(get_maps)
        .pipe(get_confidence_interval)
        .assign(region="transcript")
    )
    region = (
        get_ps(df, grouping=["csq", "region"])
        .pipe(pred_ps, poly1d=p)
        .pipe(get_maps)
        .pipe(get_confidence_interval)
    )

    # Combine the transcript- and region-level results
    maps = pd.concat([transcript, region]).sort_values(["csq", "region"])

    # Write to output
    maps.to_csv(_FILE_OUT, sep="\t", index=False)

    return maps


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
