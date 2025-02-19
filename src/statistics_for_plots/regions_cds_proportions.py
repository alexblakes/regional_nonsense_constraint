"""Get the proportion of the CDS occupied by each NMD region."""



import pandas as pd

import src
from src import statistics_for_plots as sp

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/final/nmd_annotations_simple.tsv.gz"
_FILE_OUT = "data/statistics/regions_cds_proportions.tsv"

logger = src.logger


def read_nmd_regions(path):
    logger.info("Reading NMD annotation.")

    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "pos", "enst", "region"],
        usecols=["region"],
        dtype="category",
    ).loc[:, "region"]


def proportion_value_counts(series):
    return series.value_counts(normalize=True).rename("proportion")


def append_full_cds(series):
    """Append the proportion of sites in the full CDS (1)."""
    cds = pd.Series(data=[1.0], index=["transcript"], name="proportion")
    return pd.concat([series, cds])


def write_out(series, path):
    series.to_csv(path, sep="\t")
    return series


def main():
    """Run as script."""

    cds_proportions = (
        read_nmd_regions(_FILE_IN)
        .pipe(proportion_value_counts)
        .pipe(append_full_cds)
        .pipe(sp.sort_index)
        .pipe(write_out, _FILE_OUT)
    )
    logger.info(f"CDS proportions:\n{cds_proportions}")
    logger.info("Done.")

    return cds_proportions


if __name__ == "__main__":
    src.add_log_handlers()
    main()
