"""Filter for ClinVar variants within MANE Select transcripts."""

import pandas as pd

import src

FILE_MANE = "data/interim/mane_select_transcript_ids.tsv"
FILE_CLINVAR_PARSED = "data/interim/clinvar_parsed.tsv"
FILE_OUT = "data/interim/clinvar_variants_selected.tsv"

logger = src.logger


def read_mane_transcripts(path=FILE_MANE):
    return pd.read_csv(path, sep="\t")


@src.log_step
def read_clinvar_data(path=FILE_CLINVAR_PARSED):
    return pd.read_csv(path, sep="\t")


@src.log_step
def merge_with_mane_transcript_ids(left, right):
    df = left.merge(right, how="inner", validate="many_to_one")

    logger.info(
        f"Duplicated by chr/pos/ref/alt: "
        f"{df.duplicated('chr pos ref alt'.split()).sum()}"
    )

    return df


def main():
    """Run as script."""

    mane_transcript_ids = read_mane_transcripts()

    return (
        read_clinvar_data()
        .pipe(merge_with_mane_transcript_ids, mane_transcript_ids)
        .pipe(src.write_out, FILE_OUT)
    )


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
