"""Find the number of PTVs per region, relative to the total size of the region."""



import pandas as pd

import src
from src import statistics_for_plots as sp
from src.statistics_for_plots import clinvar_proportion_vus_by_region as vus
from src.visualisation import regions_cds_proportions as cds

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/interim/clinvar_variants_vep_tidy.tsv"
_FILE_CDS_PROPORTION = "data/statistics/regions_cds_proportions.tsv"
_FILE_OUT = "data/statistics/clinvar_ptv_ascertainment.tsv"

logger = src.logger


def count_ptvs_per_region(df):
    return df.groupby("region").csq.count().rename("n_ptvs")


def count_ptvs_per_cds(df):
    return pd.Series(data=[df.csq.count()], index=["transcript"], name="n_ptvs")


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
        .pipe(sp.sort_index)
        .pipe(normalise_to_cds_footprint)
        .pipe(normalise_to_ptvs_in_cds)
        .pipe(write_out, _FILE_OUT)
    )

    logger.info(f"VUS proportions by region:\n{ptv_ascertainment}")
    logger.info("Done.")

    return ptv_ascertainment

if __name__ == "__main__":
    src.add_log_handlers()
    main()
