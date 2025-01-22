"""Boilerplate code for most modules."""

import pandas as pd
import src

FILE_IN = "data/raw/MANE.GRCh38.v0.95.summary.txt.gz"
FILE_OUT = "data/interim/mane_select_transcript_ids.tsv"
COLUMNS = {
    "symbol": "symbol",
    "RefSeq_nuc": "refseq_transcript",
    "Ensembl_nuc": "enst",
    "MANE_status": "mane_status",
}
logger = src.logger


@src.log_step
def read_data(path=FILE_IN):
    return pd.read_csv(path, sep="\t", usecols=COLUMNS).rename(columns=COLUMNS)


@src.log_step
def filter_mane_select(df):
    return df.loc[df.mane_status == "MANE Select"].drop("mane_status", axis=1)


def drop_version(df, column: str):
    df = df.assign(**{column: lambda x: x[column].str.extract(r"^(.*)\.")})
    logger.info(f"Duplicates in {column} column: {df.duplicated(column).sum()}")
    return df


def main():
    """Run as script."""

    return (
        read_data()
        .pipe(filter_mane_select)
        .pipe(drop_version, "enst")
        .pipe(drop_version, "refseq_transcript")
        .pipe(src.write_out, FILE_OUT)
    )


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
