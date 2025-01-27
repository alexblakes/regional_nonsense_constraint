"""
Find variants which are 'unexpected', e.g. pathogenic variants in unconstrained regions.

The variants are limited to those in dominant disease genes, for simplicity.
"""

import pandas as pd

import src

FILE_IN = "data/interim/clinvar_variants_in_dominant_genes.tsv.gz"
FILE_OUT = "data/interim/clinvar_unexpected_variants_by_acmg_and_constraint.tsv"
FILE_OUT_SAMPLE = (
    "data/interim/clinvar_unexpected_variants_by_acmg_and_constraint_sample.tsv"
)
logger = src.logger


@src.log_step
def read_data(path=FILE_IN):
    return pd.read_csv(path, sep="\t")


@src.log_step
def filter_ptvs(df):
    nonsense = df["csq"] == "stop_gained"
    frameshift = df["csq"] == "frameshift_variant"

    df = df[nonsense | frameshift]

    logger.info(
        f"Data shape: {df.shape}\n\n"
        f"Consequence value counts:\n{df.csq.value_counts()}"
    )

    return df


@src.log_step
def get_unexpected_variants(df):
    pathogenic = df["acmg"] == "P/LP"
    benign = df["acmg"] == "B/LB"
    constrained = df["constraint"] == "constrained"
    unconstrained = df["constraint"] == "unconstrained"

    return df.loc[(pathogenic & unconstrained) | (benign & constrained)]


@src.log_step
def reorder_columns(df):
    return df[
        "chr pos ref alt ensg enst symbol region csq acmg constraint review".split()
    ]


def random_sample(df):
    return df.groupby(["acmg", "constraint"]).sample(10)


def main():
    """Run as script."""
    return (
        read_data()
        .pipe(filter_ptvs)
        .pipe(get_unexpected_variants)
        .pipe(reorder_columns)
        .pipe(src.write_out, FILE_OUT)
        .pipe(random_sample)
        .pipe(src.write_out, FILE_OUT_SAMPLE)
    )


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
