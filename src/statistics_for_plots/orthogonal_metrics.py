"""Prepare phyloP, AlphaMissense, and pext statistics for plotting."""

from collections import defaultdict
import logging
from pathlib import Path

import pandas as pd

from src import constraint, statistics_for_plots as sp
import src
from src import constants as C

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/interim/cds_sites_phylop_pext_missense.tsv"
_FILE_OUT = "data/statistics/orthogonal_metrics.tsv.gz"
_DTYPES = defaultdict(lambda: "float16", region="category", constraint="category")

logger = logging.getLogger(__name__)


def read_data(path: str) -> pd.DataFrame:
    logger.info("Reading data.")

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["region", "constraint", "phylop", "alpha_mis", "pext"],
        dtype=_DTYPES,
    )

    logger.info(f"Valid scores:\n{df.notna().sum()}")

    return df


def tidy_constraint_values(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.query("constraint != 'indeterminate'")
        .dropna(subset="constraint")
        .assign(constraint=lambda x: x.constraint.cat.remove_unused_categories())
        .assign(
            constraint=lambda x: x.constraint.cat.rename_categories(
                {"constrained": "Constrained", "unconstrained": "Unconstrained"}
            )
        )
    )


def main():
    """Run as script."""

    annotations = (
        read_data(_FILE_IN)
        .pipe(tidy_constraint_values)
        .pipe(sp.sort_column, "region")
        .to_csv(_FILE_OUT, sep="\t", index=False)
    )

    return annotations


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
