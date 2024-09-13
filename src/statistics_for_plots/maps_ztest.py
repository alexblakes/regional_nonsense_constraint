"""Get plot statistics for MAPS."""

import logging
import itertools

import pandas as pd
from statsmodels.stats import proportion

import src
from src import constants as C
from src import statistics_for_plots as sp

_FILE_IN = "data/statistics/maps.tsv"
_FILE_OUT = "data/statistics/maps_ztest.tsv"

logger = logging.getLogger(__name__)


def read_data(path=_FILE_IN):
    return pd.read_csv(
        path,
        sep="\t",
        index_col="csq",
        usecols=["csq", "n_variants", "n_singletons", "pred_ps", "maps"],
    )


def adjust_data(df):
    return df.assign(
        ps_adj=lambda x: x["maps"] + x.loc["Synonymous (Full CDS)", "pred_ps"],
        n_singletons_adj=lambda x: x.ps_adj * x.n_variants,
    )


def z_test_one_pair(df, pair):
    df = df.loc[pair, :]
    result = proportion.proportions_ztest(df["n_singletons_adj"], df["n_variants"])
    return pd.Series(
        data=[pair[0], pair[1], result[0], result[1]], index=["left", "right", "z", "p"]
    )


def z_test_all_pairs(df):
    pairs = itertools.combinations(df.index, 2)
    results = [z_test_one_pair(df, p) for p in pairs]
    results = pd.concat(results, axis=1).T

    logger.info(f"Pairwise Z tests:\n{results}")
    logger.info(
        f"Number of tests: {len(results)}\n"
        f"Bonferroni p (0.05 / {len(results)}): {0.05 / len(results)}"
    )

    return results


def write_out(df, path):
    df.to_csv(path, sep="\t")
    return df


def main():
    """Run as script."""

    return (
        read_data().pipe(adjust_data).pipe(z_test_all_pairs).pipe(write_out, _FILE_OUT)
    )


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
