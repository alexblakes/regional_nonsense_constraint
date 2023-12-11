""" 
Get gene IDs and transcript IDs for canonical protein-coding transcripts in GENCODE v39. 
"""

# Imports
from pathlib import Path

import gtfparse  # * read_gtf makes a call to logging.basicConfig() which overwrites my logging config.
import pandas as pd

from src import constants as C
from src import setup_logger
from src.data import canonical_cds as ccds

# Logging
logger = setup_logger(Path(__file__).stem)


def get_gene_ids(df):
    """Get gene and transcript IDs, without version numbers."""

    drop_version = lambda x: x.str.split(".").str[0]

    df = (
        df[["gene_id", "transcript_id", "gene_name"]]
        .drop_duplicates()
        .apply(drop_version, axis=1)
        .drop_duplicates()
    )

    # Logging
    logger.info(f"Unique gene ids: {df.gene_id.nunique()}")
    logger.info(f"Unique transcript ids: {df.transcript_id.nunique()}")
    logger.info(f"Unique gene names: {df.gene_name.nunique()}")
    logger.info(f"Duplicated gene ids: {df.duplicated('gene_id').sum()}")
    logger.info(f"Duplicated gene ids: {df.duplicated('transcript_id').sum()}")
    logger.info(f"Duplicated gene names:\n{df[df.duplicated('gene_name', keep=False)]}")

    return df


def main():
    """Run the script."""

    df = (
        gtfparse.read_gtf(C.GENCODE_GTF, features="CDS")
        .pipe(ccds.get_canonical)
        .pipe(get_gene_ids)
    )

    # Write to output.
    logger.info("Writing gene IDs to output.")
    df.to_csv(C.CANONICAL_CDS_GENE_IDS, sep="\t", index=False)

    logger.info("Writing transcript IDs to output.")
    df.transcript_id.to_csv(
        C.CANONICAL_CDS_TRANSCRIPT_IDS, sep="\t", index=False, header=None
    )


if __name__ == "__main__":
    main()
