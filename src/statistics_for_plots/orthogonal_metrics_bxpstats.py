"""Get boxplot statistics for orthogonal metrics."""

import logging
from collections import defaultdict
import pandas as pd

import src
from src import constants as C
from src import statistics_for_plots as sfp

_CADD = "data/statistics/orthogonal_metrics_cadd.tsv.gz"
_METRICS = "data/statistics/orthogonal_metrics.tsv.gz"
_FILE_OUT = "data/statistics/orthogonal_metrics_bxp_stats.tsv"
_DTYPES = defaultdict(lambda: "float32", region="category", constraint="category")
_METRIC_NAMES = dict(
    phylop="phyloP",
    alpha_mis="AlphaMissense",
    cadd_phred="CADD Phred",
    pext="pext",
)

logger = logging.getLogger(__name__)


def read_data(path=_CADD):
    return (
        pd.read_csv(
            path,
            sep="\t",
            dtype=_DTYPES,
            # nrows=100000
        )
        .melt(id_vars=["region", "constraint"], var_name="metric", value_name="score")
        .astype({"metric": "category"})
        .pipe(sfp.sort_column, "metric", _METRIC_NAMES)
    )


def parse_cadd(path=_CADD):
    cadd = read_data(_CADD)
    cadd_full_cds = cadd.copy().assign(region="Full CDS")
    return pd.concat([cadd, cadd_full_cds]).astype({"region": "category"})


def write_out(df, path=_FILE_OUT):
    df.to_csv(path, sep="\t")

    return df


def main():
    """Run as script."""

    return (
        pd.concat([parse_cadd(_CADD), read_data(_METRICS)])
        .pipe(sfp.sort_column, labels=C.REGION_LABELS)
        .pipe(sfp.sort_column, "constraint", ["Unconstrained", "Constrained"])
        .groupby(["metric", "region", "constraint"])
        .agg(
            med=("score", lambda x: x.quantile([0.5])),
            q1=("score", lambda x: x.quantile([0.25])),
            q3=("score", lambda x: x.quantile([0.75])),
            whislo=("score", lambda x: x.quantile([0.05])),
            whishi=("score", lambda x: x.quantile([0.95])),
            mean=("score", "mean"),
            count=("score", "count"),
            std=("score","std")
        )
    ).pipe(write_out)


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
