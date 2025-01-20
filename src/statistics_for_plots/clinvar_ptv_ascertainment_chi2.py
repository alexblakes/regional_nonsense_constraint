"""Find the number of PTVs per region, relative to the total size of the region."""


import itertools

import pandas as pd
from scipy import stats

import src
from src import statistics_for_plots as sp
from src.statistics_for_plots import clinvar_proportion_vus_by_region as vus
from src.visualisation import regions_cds_proportions as cds

_FILE_IN = "data/interim/clinvar_variants_vep_tidy.tsv"
_FILE_CDS_PROPORTION = "data/statistics/regions_cds_proportions.tsv"
_FILE_OUT = "data/statistics/clinvar_ptv_ascertainment_chi2.tsv"

logger = src.logger


def count_ptvs_per_region(df):
    return df.groupby("region").csq.count().rename("n_ptvs")


def append_cds_footprint(series):
    cds_footprint = cds.read_cds_proportions(_FILE_CDS_PROPORTION).rename(
        "proportion_cds"
    )
    return pd.concat([series, cds_footprint], axis=1, join="inner")


def get_n_exp(df):
    return df.assign(n_exp=lambda x: x["n_ptvs"].sum() * x["proportion_cds"])


def tidy_for_chi2(df):
    return df.drop("proportion_cds", axis=1).rename(columns={"n_ptvs": "n_obs"})


def log_chi2_gof(df):
    logger.info(f"PTV ascertainment by region:\n{df}")
    logger.info(
        f"Chi squared goodness of fit:\n{stats.chisquare(df['n_obs'], df['n_exp'])}"
    )
    return df


def chi2_stat(df):
    o = df["n_obs"]
    e = df["n_exp"]
    return ((o - e) ** 2) / e


def log_chi2_gof_contributions(df):
    contribs = df.assign(X2=lambda x: chi2_stat(x))
    logger.info(f"Chi squared goodness of fit contributions:\n{contribs}")
    return df


def scale_n_exp(df):
    obs_totals = df["n_obs"].sum()
    exp_totals = df["n_exp"].sum()
    oe_ratio = obs_totals / exp_totals
    scaled_n_exp = df["n_exp"] * oe_ratio

    return df.assign(n_exp=scaled_n_exp)


def chi_square_one_pair(df, pair_tuple):
    df = df.loc[pair_tuple, :].pipe(scale_n_exp)
    result = stats.chisquare(df["n_obs"], df["n_exp"])

    return pd.Series(
        index=["left", "right", "X2", "p"],
        data=[pair_tuple[0], pair_tuple[1], result.statistic, result.pvalue],
    )


def chi_square_all_pairs(df):
    pairs = list(itertools.combinations(df.index, 2))
    pairwise_results = [chi_square_one_pair(df, p) for p in pairs]
    pairwise_results = pd.concat(pairwise_results, axis=1).T

    logger.info(f"Pairwise Chi square goodness of fit results:\n{pairwise_results}")

    return pairwise_results


def write_out(df, path):
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""

    ptv_ascertainment = (
        vus.read_clinvar_variants(_FILE_IN)
        .pipe(vus.filter_for_ptvs)
        .pipe(count_ptvs_per_region)
        .pipe(sp.sort_index)
        .pipe(append_cds_footprint)
        .drop("Full CDS")
        .pipe(get_n_exp)
        .pipe(tidy_for_chi2)
        .pipe(log_chi2_gof)
        .pipe(log_chi2_gof_contributions)
        .pipe(chi_square_all_pairs)
        .pipe(write_out, _FILE_OUT)
    )

    logger.info("Done.")

    return ptv_ascertainment


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
