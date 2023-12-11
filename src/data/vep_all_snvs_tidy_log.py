"""Log the number of variants annotated by VEP."""

# Imports
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C

# Logging
logger = setup_logger(Path(__file__).stem)


# Module constants


# Functions
def count_vep_vars(path):
    """Count variants in the raw VEP output."""

    logger.info(f"Counting variants with VEP annotation.")

    n = len(
        pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            usecols=[1],
            dtype="category",
        )
    )
    return n


def get_vep_tidy(path):
    """Read the tidied VEP data."""

    logger.info(f"Reading tidied VEP data.")

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "pos", "ref", "alt", "csq", "enst"],
        dtype={
            "chr": "category",
            "pos": "int32",
            "ref": "category",
            "alt": "category",
            "csq": "category",
            "enst": "category",
        },
    )

    return df


def main():
    """Run the script."""

    # Read data
    vep = get_vep_tidy(C.VEP_ALL_SNVS_TIDY)
    gene_ids = pd.read_csv(C.CANONICAL_CDS_GENE_IDS, sep="\t", dtype="category")

    # Log key results
    logger.info(f"Number of SNVs after tidying: {len(vep)}")
    logger.info(f"Duplicated SNVs: {vep.duplicated().sum()}")
    logger.info("NB all possible SNVs in the CDS include start codon variants, which are not recorded in the VEP output.")
    logger.info(
        f"Number of canonical transcripts in VEP annotation: {vep.enst.nunique()}"
    )
    logger.info(
        f"Variant counts in canonical transcripts:\n{vep.csq.value_counts()}"
    )


if __name__ == "__main__":
    main()
