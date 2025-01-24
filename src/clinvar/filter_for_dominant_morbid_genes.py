"""Filter ClinVar variants for those in OMIM morbid dominant genes."""

import pandas as pd
import src

FILE_CLINVAR = "data/interim/clinvar_variants_constraint.tsv.gz"
FILE_OMIM = "data/interim/genemap2_simple.tsv"
FILE_OUT = "data/interim/clinvar_variants_constraint_in_morbid_dominant_genes.tsv.gz"
logger = src.logger


@src.log_step(display_args=True)
def read_data(path):
    return pd.read_csv(path, sep="\t")


@src.log_step
def filter_ad_morbid_genes(df):
    ad_mask = df["inheritance"] == "Autosomal dominant"
    return df.loc[ad_mask, ["ensg"]].drop_duplicates()


@src.log_step
def merge(left, right):
    return left.merge(right, how="inner", validate="many_to_one")


def main():
    """Run as script."""

    omim = read_data(path=FILE_OMIM).pipe(filter_ad_morbid_genes)
    
    clinvar = (
        read_data(path=FILE_CLINVAR).pipe(merge, omim).pipe(src.write_out, FILE_OUT)
    )

    return clinvar


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
