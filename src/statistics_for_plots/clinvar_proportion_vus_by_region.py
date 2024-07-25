"""Find the proportion of VUS PTVs per region."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import statistics_for_plots as sp

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/interim/clinvar_variants_vep_tidy.tsv"
_FILE_OUT = "data/statistics/clinvar_vus_by_region.tsv"

logger = logging.getLogger(__name__)


def read_clinvar_variants(path):
    logger.info("Reading ClinVar data")
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "pos", "enst", "ref", "alt", "csq", "region", "acmg"],
        usecols=["csq", "region", "acmg"],
        dtype="category",
        na_values=".",
    )


def filter_for_ptvs(df):
    return df.query("csq in ['stop_gained', 'frameshift']")


def get_proportion_vus_per_cds(df):
    return (
        df.acmg.value_counts(normalize=True)
        .rename("proportion_vus")
        .loc[["VUS"]]
        .set_axis(["transcript"])
    )


def get_proportion_vus_per_region(df):
    return (
        df.groupby("region")
        .acmg.value_counts(normalize=True)
        .rename("proportion_vus")
        .xs("VUS", level="acmg")
    )


def concat_proportion_vus_in_cds_and_regions(df):
    return pd.concat(
        [get_proportion_vus_per_cds(df), get_proportion_vus_per_region(df)]
    )


def write_out(series, path):
    series.to_csv(path, sep="\t")
    return series


def main():
    """Run as script."""

    vus_proportions = (
        read_clinvar_variants(_FILE_IN)
        .pipe(filter_for_ptvs)
        .pipe(concat_proportion_vus_in_cds_and_regions)
        .pipe(sp.sort_index)
        .pipe(write_out, _FILE_OUT)
    )

    logger.info(f"VUS proportions by region:\n{vus_proportions}")
    logger.info("Done.")

    return


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
