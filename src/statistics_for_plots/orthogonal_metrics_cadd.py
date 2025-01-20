"""Prepare cadd statistics for plotting."""

from collections import defaultdict


import pandas as pd

from src import statistics_for_plots as sp
from src.statistics_for_plots import orthogonal_metrics
import src

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/interim/cadd_scores_coding_annotated.tsv"
_FILE_OUT = "data/statistics/orthogonal_metrics_cadd.tsv.gz"
_DTYPES = defaultdict(lambda: "float16", region="category", constraint="category")

logger = src.logger


def read_data(path: str) -> pd.DataFrame:
    logger.info("Reading data.")

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["region", "constraint", "cadd_phred"],
        dtype=_DTYPES,
    )

    logger.info(f"Valid scores:\n{df.notna().sum()}")

    return df


def main():
    """Run as script."""

    cadd = (
        read_data(_FILE_IN)
        .pipe(orthogonal_metrics.tidy_constraint_values)
        .pipe(sp.sort_column, "region")
        .to_csv(_FILE_OUT, sep="\t", index=False)
    )

    return cadd


if __name__ == "__main__":
    src.add_log_handlers()
    main()
