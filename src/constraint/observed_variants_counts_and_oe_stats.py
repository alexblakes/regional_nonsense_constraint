"""Get variant counts and O/E statistics by variant consequence and region.

Constraint statistics are only calculated using variants above a coverage level.
"""

import argparse
import logging
from pathlib import Path

from fast_poibin import PoiBin
import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

import src
from src import constants as C

_FILE_IN = "data/interim/cds_all_possible_snvs_annotated.tsv.gz"
_FILE_OUT = "data/interim/oe_stats_regions_cov_"  # Suffix appended when written
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
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
_MU = "mu_roulette_scaled"

logger = logging.getLogger(__name__)


def get_variant_annotations(path):
    """Read variant annotation data."""

    logger.info("Reading variant annotation data.")

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=_NAMES_DICT,
        comment="#",
        header=None,
        names=_NAMES_DICT.keys(),
        na_values=".",
        # nrows=100000,  #! Testing
    )

    return df


def mark_observed_variants(df):
    """Assign a new boolean column showing observed variants."""

    df = df.fillna({"ac": 0})
    df["obs"] = np.where(df["ac"] >= 1, True, False)

    return df


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


def per_row_binom_ci(row):
    try:
        result = stats.binomtest(
            row["n_obs"], row["n_pos"], row["prop_exp"], alternative="less"
        )

        # Return the upper bound of the 95% CI around the O/E value
        return (result.proportion_ci().high * row["n_pos"]) / row["n_exp"]

    except:
        pass


def per_row_poisson_binom_p(row):
    try:
        poibin = PoiBin(row["mu_list"])

        # Return a one-sided P value (H0: n_obs >= n_exp)
        return poibin.cdf[row["n_obs"]]

    except:
        pass


def agg_stats(df, grouping_columns, mu=_MU):
    """Get variant counts and O/E statistics for a given grouping."""

    tqdm.pandas(desc=f"Aggregate statistics, grouping by {grouping_columns}")

    df = (
        df.groupby(grouping_columns, observed=True)
        .agg(
            n_obs=("obs", "sum"),
            n_pos=("pos", "count"),
            n_exp=(mu, "sum"),
            mu_list=(mu, lambda x: list(x)),
        )
        .assign(
            oe=lambda x: x.n_obs / x.n_exp,
            prop_obs=lambda x: x.n_obs / x.n_pos,
            prop_exp=lambda x: x.n_exp / x.n_pos,  # Equals mu.mean()
            oe_ci_hi=lambda x: x.progress_apply(per_row_binom_ci, axis=1),
            p=lambda x: x.apply(per_row_poisson_binom_p, axis=1),
        )
        .drop("mu_list", axis=1)
        .reset_index()
    )

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


def parse_args():
    """Parse command line arguments."""

    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument(
        "-c",
        "--coverage",
        type=int,
        default=[20],
        nargs="*",
        help="Minimum coverage of sites to include. Accepts multiple integer values.",
    )

    return parser.parse_args()


def main():
    """Run as script."""

    df = get_variant_annotations(_FILE_IN).pipe(mark_observed_variants)

    # Filter for coverage
    coverage = parse_args().coverage
    for cov in coverage:
        df_min_coverage = filter_covered_sites(df, cov)

        # Get variants by region for constraint calculations.
        # ? Should we limit to rare variants only?
        transcript = (
            agg_stats(df_min_coverage, ["enst", "csq"])
            .assign(region="transcript")
            .pipe(reorder_data)
        )

        regions = agg_stats(df_min_coverage, ["enst", "csq", "region"]).pipe(
            reorder_data
        )

        combined = pd.concat([regions, transcript])

        # Logging
        logger.info(f"Counts grouped by enst/csq: {len(transcript)}")
        logger.info(f"Counts grouped by nmd/enst/csq: {len(regions)}")
        logger.info(f"Number of entries: {len(combined)}")
        logger.info(f"NaN values:\n{combined.isna().sum()}")
        logger.info(f"Unique transcripts: {combined.enst.nunique()}")
        logger.info(
            f"Region value counts:\n{combined.groupby('region').csq.value_counts()}"
        )

        # Write to output
        logger.info("Writing to output.\n")
        combined.to_csv(f"{_FILE_OUT}{str(cov)}.tsv", sep="\t", index=False)

    return combined  #! Testing


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
