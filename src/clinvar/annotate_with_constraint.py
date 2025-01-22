"""Add constraint annotation to ClinVar variants."""

import pandas as pd

import src

CLINVAR = "data/interim/clinvar_variants_vep.tsv.gz"
CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
FILE_OUT = "data/interim/clinvar_variants_constraint.tsv.gz"

logger = src.logger


@src.log_step
def read_clinvar(path=CLINVAR):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names="chr pos ref alt symbol enst csq region acmg review".split(),
        na_values=".",
    )


@src.log_step
def read_constraint(path):
    return pd.read_csv(path, sep="\t", usecols=["enst", "region", "constraint"])


def merge(left, right):
    df = left.merge(right, how="inner", validate="many_to_one")
    logger.info(f"Data shape: {df.shape}\n\n" f"NaN values:\n{df.isna().sum()}")
    return df


def main():
    """Run as script."""
    constraint = read_constraint(CONSTRAINT)
    return read_clinvar(CLINVAR).pipe(merge, constraint).pipe(src.write_out, FILE_OUT)


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
