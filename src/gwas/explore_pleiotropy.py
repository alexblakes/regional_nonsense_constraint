"""Exploring pleiotropy in GWAS catalogue."""

import pandas as pd

import src

FILE_IN = "data/interim/gwas_catalog_tidy.tsv.gz"
FILE_OUT = "data/statistics/gwas_pleiotropy_by_csq.tsv"
logger = src.logger


def main():
    """Run as script."""
    df = src.read_data(FILE_IN, low_memory=False)

    variant_counts = (
        df.drop_duplicates(["chrom", "pos", "alt"])
        .csq.value_counts()
        .rename("unique_variants")
    )

    trait_counts = df.groupby("csq").agg(unique_traits=("trait", "nunique"))

    return (
        pd.concat([variant_counts, trait_counts], axis=1)
        .assign(traits_per_variant=lambda x: x["unique_traits"] / x["unique_variants"])
        .sort_values("traits_per_variant", ascending=False)
        .pipe(src.write_out, FILE_OUT)
    )


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
