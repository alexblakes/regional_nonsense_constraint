"""Filter ClinVar variants for those in OMIM morbid dominant genes."""

import pandas as pd
import src

FILE_CLINVAR = "data/interim/clinvar_variants_constraint.tsv.gz"
FILE_OMIM = "data/interim/genemap2_simple.tsv"
logger = src.logger


@src.log_step(display_args=True)
def read_data(path):
    return pd.read_csv(path, sep="\t")


def main():
    """Run as script."""

    clinvar = read_data(path=FILE_CLINVAR)
    omim = read_data(path=FILE_OMIM)

    return clinvar, omim


if __name__ == "__main__":
    src.add_log_handlers()
    clinvar, omim = main()
