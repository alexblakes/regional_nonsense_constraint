"""Get the proportion of the CDS occupied by each NMD region."""

import logging
from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/final/nmd_annotations_simple.tsv.gz"
_FILE_OUT = "data/statistics/regions_cds_proportions.tsv"
_CATEGORIES = ["full_cds", "nmd_target", "start_proximal", "long_exon", "distal_nmd"]
_LABELS = ["Full CDS", "NMD target", "Start proximal", "Long exon", "Distal"]

logger = logging.getLogger(__name__)


def read_nmd_regions(path):
    logger.info("Reading NMD annotation.")

    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "pos", "enst", "region"],
        nrows=100000,
        usecols=["region"],
        dtype="category",
    ).loc[:, "region"]


def proportion_value_counts(series):
    return series.value_counts(normalize=True).rename("proportion")


def append_full_cds(series):
    """Append the proportion of sites in the full CDS (1)."""
    cds = pd.Series(data=[1.0], index=["full_cds"], name="proportion")
    return pd.concat([series, cds])


def sort_index(series, categories, labels, **kwargs):
    """Create an ordered categorical index for a series."""
    _ix = pd.CategoricalIndex(categories, **kwargs)
    return series.reindex(_ix).rename(dict(zip(categories, labels)))


def write_out(series, path):
    series.to_csv(path, sep="\t")
    return series


def main():
    """Run as script."""

    cds_proportions = (
        read_nmd_regions(_FILE_IN)
        .pipe(proportion_value_counts)
        .pipe(append_full_cds)
        .pipe(sort_index, _CATEGORIES, _LABELS, ordered=True, name="region")
        .pipe(write_out, _FILE_OUT)
    )
    logger.info(f"CDS proportions:\n{cds_proportions}")
    logger.info("Done.")

    return cds_proportions


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
