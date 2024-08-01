"""Get variant counts and O/E statistics by variant consequence and region.

Constraint statistics are only calculated using variants above a coverage level.
"""

import logging
from pathlib import Path

from fast_poibin import PoiBin
import numpy as np
import pandas as pd
from scipy import stats
from sklego import pandas_utils as pu
from tqdm import tqdm

import src

_FILE_IN = "data/interim/cds_all_possible_snvs_annotated.tsv.gz"
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_MU = "mu_roulette_scaled"  # The mutation rate estimate used in constraint calculations
_MIN_COVERAGE = 20
_FILE_OUT = f"data/interim/oe_stats_regions_cov_{_MIN_COVERAGE}.tsv"
_NAMES_DICT = {
    "chr": "category",
    "pos": "int32",
    "ref": "category",
    "alt": "category",
    "enst": "category",
    "csq": "category",
    "region": "category",
    "median_coverage": "Int16",
    "ac": "Int32",
    "an": "Int32",
    "af": "Float32",
    "pent": "category",
    "mu_roulette": "float32",
    "mu_roulette_final": "float32",
    "mu_roulette_scaled": "float32",
    "mu_gnomad": "float32",
}

logger = logging.getLogger(__name__)


def read_data(path):
    """Read variant annotation data."""

    logger.info("Reading variant annotation data.")

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=_NAMES_DICT,
        comment="#",
        header=None,
        names=_NAMES_DICT.keys(), # type: ignore
        na_values=".",
    )

    return df


def mark_observed_variants(df):
    """Assign a new boolean column showing observed variants."""

    return df.fillna({"ac": 0}).assign(
        obs=lambda x: np.where(x["ac"] >= 1, True, False)
    )


@pu.log_step(print_fn=lambda x: logger.info(x), shape_delta=True)
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
            f"Qualifying variants by observed / consequence:\n{df.groupby('obs').csq.value_counts()}"
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


def drop_mu_nans(grouped_series):
    """Drop invalid mutation rate annotations (NaNs) from the series.
    
    Some variants do not have a mutability score. For regions harbouring these sites,
    the number of expected variants will be under-estimated. The constraint results will
    be more conservative in these regions.
    """
    return [x for x in grouped_series if 0 <= x <= 1]


def per_row_binom_test(row):
    return stats.binomtest(
        row["n_obs"], row["n_pos"], row["prop_exp"], alternative="less"
    )


def per_row_binom_oe_ci_hi(row):
    """Find the upper bound of the 95% CI around the O/E value."""
    if row["n_exp"] == 0: # binomtest not applicable.
        return np.nan
    
    return (per_row_binom_test(row).proportion_ci().high * row["n_pos"]) / row["n_exp"]


def per_row_poisson_binom_p(row):

    # Rarely, there seem to be more observed variants than possible variant sites. This
    # is because sites without a mutability score have been dropped. We can't run a 
    # Poissin-Binomial test in this case.
    if row["n_obs"] > len(row["mu_list"]):
        return np.nan
    
    return PoiBin(row["mu_list"]).cdf[row["n_obs"]]


def agg_stats(df, grouping_columns, mu=_MU):
    """Get variant counts and O/E statistics for a given grouping."""

    tqdm.pandas(desc=f"Aggregate statistics, grouping by {grouping_columns}")

    df = (
        df.groupby(grouping_columns, observed=True)
        .agg(
            n_obs=("obs", "sum"),
            n_pos=("pos", "count"),
            n_exp=(mu, "sum"),
            mu_list=(mu, drop_mu_nans),
        )
        .assign(
            oe=lambda x: x.n_obs / x.n_exp,
            prop_obs=lambda x: x.n_obs / x.n_pos,
            prop_exp=lambda x: x.n_exp / x.n_pos,  # Equals mu.mean()
            oe_ci_hi=lambda x: x.progress_apply(per_row_binom_oe_ci_hi, axis=1),
            p=lambda x: x.apply(per_row_poisson_binom_p, axis=1),
        )
        .drop("mu_list", axis=1)
        .reset_index()
    )

    return df


def get_stats_per_region(df):
    df = agg_stats(df, ["enst", "csq", "region"])
    logger.info(f"Counts grouped by enst/csq/region: {len(df)}")
    return df


def get_stats_per_transcript(df):
    df = agg_stats(df, ["enst", "csq"]).assign(region="transcript")
    logger.info(f"Counts grouped by enst/csq: {len(df)}")
    return df


def combine_stats(df):
    df = pd.concat([get_stats_per_region(df), get_stats_per_transcript(df)])

    logger.info(f"Number of entries: {len(df)}")
    logger.info(f"NaN values:\n{df.isna().sum()}")
    logger.info(f"Unique transcripts: {df.enst.nunique()}")
    logger.info(f"Region value counts:\n{df.groupby('region').csq.value_counts()}")

    return df


def reorder_data(df):
    return df[
        [
            "enst",
            "region",
            "csq",
            "n_pos",
            "n_obs",
            "n_exp",
            "prop_obs",
            "prop_exp",
            "oe",
            "oe_ci_hi",
            "p",
        ]
    ]


def write_out(df, path):
    logger.info("Writing to output.\n")
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""

    return (
        read_data(_FILE_IN)
        .pipe(mark_observed_variants)
        .pipe(filter_covered_sites, _MIN_COVERAGE)
        .pipe(combine_stats)
        .pipe(reorder_data)
        .pipe(write_out, f"{_FILE_OUT}")
    )


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
