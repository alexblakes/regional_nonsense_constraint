"""Conduct pairwise Z tests for proportions."""

import logging
import itertools

import pandas as pd
from statsmodels.stats import proportion

import src
from src.statistics_for_plots import clinvar_proportion_vus_by_region as cv

_FILE_OUT = "data/statistics/clinvar_vus_by_region_ztest.tsv"

logger = logging.getLogger(__name__)


def get_proportion_vus_per_region(df):
    return (
        df.groupby("region")
        .agg(
            n_variants=("acmg", "count"),
            n_vus=("acmg", lambda x: sum([y == "VUS" for y in x])),
        )
        .assign(proportion_vus=lambda x: x.n_vus / x.n_variants)
    )


def z_test_one_pair(df, pair):
    df = df.loc[pair, :]
    result = proportion.proportions_ztest(df["n_vus"], df["n_variants"])

    return pd.Series(data=[pair[0], pair[1], result[0], result[1]], index=["left","right","z","p"])


def z_test_all_pairs(df):
    pairs = itertools.combinations(df.index, 2)
    results = [z_test_one_pair(df, p) for p in pairs]
    results =  pd.concat(results, axis=1).T

    logger.info(f"Pairwise Z test of proportions:\n{results}")

    return results


def write_out(df, path):
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""

    vus_proportions = (
        cv.read_clinvar_variants()
        .pipe(cv.filter_for_ptvs)
        .pipe(get_proportion_vus_per_region)
        .pipe(z_test_all_pairs)
        .pipe(write_out, _FILE_OUT)
    )

    return vus_proportions


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
