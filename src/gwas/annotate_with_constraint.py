"""Annotate GWAS PTVs with constraint."""

import pandas as pd

import src
from src import constants as C
from src import statistics_for_plots as sp

FILE_GWAS_PTVS = "data/interim/gwas_ptvs_by_region.tsv"
FILE_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
FILE_OUT = "data/final/gwas_ptvs_in_constrained_regions.tsv"
FILE_STATS = "data/statistics/gwas_ptvs_in_constrained_regions.tsv"

logger = src.logger


@src.log_step
def read_gwas_ptvs(path=FILE_GWAS_PTVS):
    return src.read_data(
        path,
        header=None,
        names=["chrom", "start", "pos", "alt", "enst", "region", "csq"],
    )


@src.log_step
def read_constraint(path=FILE_CONSTRAINT):
    return src.read_data(path, usecols=["enst", "region", "constraint"])


def merge(left, right):
    df = pd.merge(left, right, how="inner", validate="many_to_one")
    logger.info(f"Data shape: {df.shape}\n\n" f"NaN values: {df.isna().sum()}")
    return df


@src.log_step
def drop_nans(df):
    return df.dropna()


def main():
    """Run as script."""

    constraint = read_constraint(FILE_CONSTRAINT)

    return (
        read_gwas_ptvs()
        .pipe(merge, constraint)
        .pipe(drop_nans)
        .pipe(src.write_out, FILE_OUT)
        .value_counts("constraint")
        .rename("count").reset_index()
        .pipe(sp.sort_column, "constraint", C.CONSTRAINT_LABELS)
        .pipe(src.write_out, FILE_STATS, index=True)
    )


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
