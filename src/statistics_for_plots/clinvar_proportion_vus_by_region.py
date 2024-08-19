"""Find the proportion of VUS PTVs per region."""

import logging
from pathlib import Path

import pandas as pd
from statsmodels.stats import proportion

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


def per_row_proportion_z_test(row, value, **kwargs):
    return proportion.proportions_ztest(row.n_vus, row.n_variants, value)[1]


def get_proportion_vus_per_region(df):
    return (
        df.groupby("region")
        .agg(
            n_variants=("acmg", "count"),
            n_vus=("acmg", lambda x: sum([y == "VUS" for y in x])),
        )
        .assign(proportion_vus=lambda x: x.n_vus / x.n_variants)
        .assign(ci_lo=lambda x: proportion.proportion_confint(x.n_vus, x.n_variants)[0])
        .assign(err=lambda x: x.proportion_vus - x.ci_lo)  # The error is symmetrical
        .assign(
            p=lambda x: x.apply(
                per_row_proportion_z_test,
                value=x.loc["nmd_target", "proportion_vus"],
                axis=1,
            )
        )
        .assign(bfr_p=lambda x: x.p < 0.05 / (len(x) - 1))
        .loc[:, ["proportion_vus", "err", "p", "bfr_p"]]
    )


def add_transcript_as_region(df):
    return pd.concat([df, df.assign(region="transcript")])


def write_out(series, path):
    series.to_csv(path, sep="\t")
    return series


def main():
    """Run as script."""

    vus_proportions = (
        read_clinvar_variants(_FILE_IN)
        .pipe(filter_for_ptvs)
        .pipe(add_transcript_as_region)
        .pipe(get_proportion_vus_per_region)
        .pipe(sp.sort_index)
        .pipe(write_out, _FILE_OUT)
    )

    logger.info(f"VUS proportions by region:\n{vus_proportions}")
    logger.info("Done.")

    return vus_proportions


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
