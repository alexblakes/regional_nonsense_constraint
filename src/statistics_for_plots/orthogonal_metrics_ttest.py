"""Perform Welch's T test for metrics in constrained vs unconstrained regions."""

import logging
import pandas as pd
from scipy import stats

import src

_FILE_IN = "data/statistics/orthogonal_metrics_bxp_stats.tsv"

logger = logging.getLogger(__name__)


def read_data(path=_FILE_IN):
    return pd.read_csv(path, sep="\t")


def main():

    df = read_data().set_index(["metric", "region", "constraint"])
    
    constrained = df.xs("Constrained", level="constraint")
    unconstrained = df.xs("Unconstrained", level="constraint")

    res = stats.ttest_ind_from_stats(
        constrained["mean"],
        constrained["std"],
        constrained["count"],
        unconstrained["mean"],
        unconstrained["std"],
        unconstrained["count"],
        equal_var=False,
        alternative="two-sided"
    )

    logger.info(f"Welch's T test statistics:\n{res}")

    return None


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
