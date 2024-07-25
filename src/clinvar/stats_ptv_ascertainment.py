"""Find the number of PTVs per region, relative to the total size of the region."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import utils
from src.clinvar import stats_proportion_vus_by_region as vus
from src.visualisation import plot_cds_proportions as cds

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/interim/clinvar_variants_vep_tidy.tsv"
_FILE_CDS_PROPORTION = "data/statistics/regions_cds_proportions.tsv"
_FILE_OUT = "data/statistics/clinvar_ptv_ascertainment.tsv"

logger = logging.getLogger(__name__)


def count_ptvs_per_region(df):
    return df.groupby("region").csq.count().rename("n_ptvs")


def count_ptvs_per_cds(df):
    return pd.Series(data=[df.csq.count()], index=["full_cds"], name="n_ptvs")


def concat_ptvs_per_region_and_cds(df):
    return pd.concat([count_ptvs_per_region(df), count_ptvs_per_cds(df)])


def normalise_to_ptvs_in_cds(series):
    return series / series.loc["Full CDS"]


def normalise_to_cds_footprint(series):
    cds_footprint = cds.read_cds_proportions(_FILE_CDS_PROPORTION)
    return (series / cds_footprint).rename("ptv_ascertainment")


def write_out(series, path):
    series.to_csv(path, sep="\t")
    return series


def main():
    """Run as script."""

    ptv_ascertainment = (
        vus.read_clinvar_variants(_FILE_IN)
        .pipe(vus.filter_for_ptvs)
        .pipe(concat_ptvs_per_region_and_cds)
        .pipe(utils.sort_index)
        .pipe(normalise_to_cds_footprint)
        .pipe(normalise_to_ptvs_in_cds)
        .pipe(write_out, _FILE_OUT)
    )

    logger.info(f"VUS proportions by region:\n{ptv_ascertainment}")
    logger.info("Done.")

    return ptv_ascertainment

if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
