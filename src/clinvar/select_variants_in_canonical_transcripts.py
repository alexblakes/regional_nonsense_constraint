"""Select ClinVar variants from canonical transcripts of interest."""

import pandas as pd
from sklego import pandas_utils

import src

FILE_CANONICAL_TRANSCRIPTS = "data/final/canonical_transcript_ids.txt"
FILE_CLINVAR_PARSED = "data/interim/clinvar_parsed.tsv"
FILE_BIOMART_TRANSCRIPT_IDS = (
    "data/raw/250121_biomart_ensembl_refseq_transcript_ids.txt.gz"
)

logger = src.logger


@pandas_utils.log_step_extra(lambda x: f"Length: {len(x)}", print_fn=logger.info)
def get_canonical_enst_ids(path=FILE_CANONICAL_TRANSCRIPTS):
    return pd.read_csv(path).squeeze().drop_duplicates().rename("enst")


def read_clinvar_data(path=FILE_CLINVAR_PARSED):
    return pd.read_csv(path, sep="\t")


@src.log_step
def read_biomart_transcript_ids(path=FILE_BIOMART_TRANSCRIPT_IDS):
    return pd.read_csv(path, sep="\t", header=0, names=["enst", "refseq_transcript"])


@src.log_step
def drop_nans(df):
    return df.dropna()


def drop_duplicates(df):
    df = df.drop_duplicates()
    logger.info(
        f"Data shape: {df.shape}\n\n"
        f"Duplicated RefSeq IDs: {df.duplicated('refseq_transcript').sum()}"
    )
    return df


def merge(left, right):
    return left.merge(right, how="inner", validate="many_to_one")


def main():
    """Run as script."""

    canonical_ensts = get_canonical_enst_ids()
    clinvar_parsed = read_clinvar_data()
    biomart_transcript_ids = (
        read_biomart_transcript_ids()
        .pipe(drop_nans)
        .pipe(merge, canonical_ensts)
        .pipe(drop_duplicates)
        # .pipe(pd.merge, clinvar_parsed, validate="one_to_many")
    )

    return biomart_transcript_ids


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
