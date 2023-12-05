"""Reformat pext scores to bed file for liftover."""

# Imports
from pathlib import Path

import pandas as pd

from src import constants as C
from src import setup_logger


# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_pext_data(path):
    """Read pext data to memory."""

    logger.info("Reading pext data.")

    df = (
        pd.read_csv(
            path,
            sep="\t",
            # nrows=100000,  #! Testing
        )
        .dropna()
        .reset_index(drop=True)
    )

    logger.info(f"pext annotations: {len(df)}")

    return df


def reformat_to_bed(df):
    """Reformat to bed."""

    logger.info("Reformating to bed.")

    locus = df.locus.str.split(":")

    df["chr"] = "chr" + locus.str[0]
    df["end"] = locus.str[1].astype("int32")
    df["start"] = df["end"] - 1

    df = df.rename(columns={"mean_proportion": "mean_pext"})

    df = df[["chr", "start", "end", "ensg", "mean_pext"]]

    logger.info(f"Unique chr values: {df.chr.unique()}")

    return df


def main():
    """Run as script."""

    df = read_pext_data(C.PEXT_RAW).pipe(reformat_to_bed)

    # Write to .bed output
    logger.info("Writing bed file to output.")
    df.to_csv(C.PEXT_BED_37, sep="\t", index=False, header=False)

    return df  #! Testing


if __name__ == "__main__":
    main()
