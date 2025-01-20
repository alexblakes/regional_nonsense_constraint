"""Calculate correlation coefficients for regional constraint vs LOEUF / pLI."""

import itertools


import pandas as pd
from scipy import stats

import src

FILE_IN = "data/final/regional_nonsense_constraint.tsv"

logger = src.logger


def read_data(path):
    return pd.read_csv(
        path, sep="\t", usecols=["region", "oe", "oe_ci_hi", "pli", "loeuf"]
    )


def main():
    """Run as script."""
    df = read_data(FILE_IN).query("region == 'transcript'").dropna()

    logger.info(f"Transcript count (NaNs dropped): {len(df)}")

    oe_scores = ["oe", "oe_ci_hi"]
    gnomad_scores = ["loeuf", "pli"]

    for a, b in itertools.product(oe_scores, gnomad_scores):
        variables = df[[a, b]]
        rho, p = stats.spearmanr(variables, nan_policy="omit")
        logger.info(f"{a} vs {b}; rho={rho}, p={p}")


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
