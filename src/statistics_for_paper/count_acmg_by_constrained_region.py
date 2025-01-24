"""Boilerplate code for most modules."""

import pandas as pd
from statsmodels.stats import proportion

import src

FILE_IN = "data/interim/clinvar_variants_constraint_in_morbid_dominant_genes.tsv.gz"
FILE_OUT = "data/statistics/clinvar_ptvs_in_morbid_ad_genes_acmg_by_constraint_and_csq.tsv"

logger = src.logger


@src.log_step
def read_data(path=FILE_IN):
    usecols = ["csq", "region", "acmg", "constraint"]
    df = pd.read_csv(path, sep="\t", usecols=usecols)

    logger.info(
        f"Data shape: {df.shape}\n\n" f"PTV value counts:\n{df.csq.value_counts()}"
    )

    return df


def count_acmg_and_constraint(df):
    return (
        df.groupby(["csq", "constraint"])["acmg"]
        .value_counts()
        .to_frame("count")
        .reset_index()
        .assign(
            total=lambda x: x.groupby(["csq", "constraint"])["count"].transform("sum"),
            proportion=lambda x: x["count"] / x["total"],
            ci95_lo=lambda x: proportion.proportion_confint(x["count"], x["total"])[0],
            err95=lambda x: x["proportion"] - x["ci95_lo"],
        )
        .drop("ci95_lo", axis=1)
    )


def main():
    """Run as script."""
    return read_data().pipe(count_acmg_and_constraint).pipe(src.write_out, FILE_OUT)


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
