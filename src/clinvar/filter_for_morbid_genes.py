"""Filter ClinVar variants for those in OMIM morbid dominant genes."""

import pandas as pd
import src

FILE_CLINVAR = "data/interim/clinvar_variants_constraint.tsv.gz"
FILE_OMIM = "data/interim/genemap2_simple.tsv"
FILE_OUT_DOMINANT = "data/interim/clinvar_variants_in_dominant_genes.tsv.gz"
FILE_OUT_RECESSIVE = "data/interim/clinvar_variants_in_recessive_genes.tsv.gz"
logger = src.logger


@src.log_step(display_args=True)
def read_data(path):
    return pd.read_csv(path, sep="\t")


@src.log_step(display_args=True)
def filter_by_inheritance(df, inheritance):
    ad_mask = df["inheritance"] == inheritance
    return df.loc[ad_mask, ["ensg"]].drop_duplicates()


@src.log_step
def merge(left, right):
    return left.merge(right, how="inner", validate="one_to_many")


@src.log_step
def exclude_overlapping(df, other, column):
    return df[~df[column].isin(other[column])]


def main():
    """Run as script."""

    omim = read_data(path=FILE_OMIM)
    clinvar = read_data(path=FILE_CLINVAR)

    dominant = (
        filter_by_inheritance(omim, "Autosomal dominant")
        .pipe(merge, clinvar)
        .pipe(src.write_out, FILE_OUT_DOMINANT)
    )
    recessive = (
        filter_by_inheritance(omim, "Autosomal recessive")
        .pipe(merge, clinvar)
        .pipe(exclude_overlapping, dominant, "ensg")
        .pipe(src.write_out, FILE_OUT_RECESSIVE)
    )

    return dominant, recessive


if __name__ == "__main__":
    src.add_log_handlers()
    ad, ar = main()
