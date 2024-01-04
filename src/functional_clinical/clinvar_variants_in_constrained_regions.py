"""Identify ClinVar LoF variants in constrained regions."""

# Imports
from pathlib import Path

import numpy as np
import pandas as pd

from src import constants as C
from src import setup_logger
from src.functional_clinical import merge_orthogonal_annotations as moa


# Logging
logger = setup_logger(Path(__file__).stem)


# Module constants
_PLOF = ["stop_gained", "frameshift_variant"]
_CONSTRAINT_COLS = ["enst", "region", "constraint", "pli", "loeuf"]


# Functions
def read_clinvar_variants(path):
    """Read the tidied, VEP-annotated ClinVar data."""

    logger.info(f"Reading ClinVar data.")

    df = pd.read_csv(path, sep="\t")

    logger.info(f"Tidied ClinVar variants: {len(df)}")
    logger.info(f"NaN values:\n{df.isna().sum()}")

    return df


def keep_truncating_variants(df):
    """Filter for frameshift and stop-gained variants."""

    df = df[df.csq.isin(_PLOF)]

    logger.info(f"pLoF variants: {len(df)}")
    logger.info(f"pLoF value counts:\n{df.csq.value_counts()}")

    return df


def simplify_acmg_annotations(df):
    """Simplify ACMG annotations."""

    logger.info(f"ACMG value counts:\n{df.acmg.value_counts()}")

    df = df.replace(
        {
            "Pathogenic": "P/LP",
            "Likely pathogenic": "P/LP",
            "Uncertain significance": "VUS",
            "Benign": "B/LB",
            "Likely benign": "B/LB",
        }
    )

    logger.info(f"Simplified ACMG value counts:\n{df.acmg.value_counts()}")

    return df


def main():
    """Run as script."""

    # Read datasets
    clinvar = (
        read_clinvar_variants(C.CLINVAR_VEP_TIDY)
        .pipe(keep_truncating_variants)
        .pipe(simplify_acmg_annotations)
    )

    nmd = moa.read_nmd_annotations(C.NMD_ANNOTATIONS)

    constraint = pd.read_csv(
        C.REGIONAL_NONSENSE_CONSTRAINT, sep="\t", usecols=_CONSTRAINT_COLS
    )

    # Annotate variants with their NMD regions
    df = clinvar.merge(nmd, how="inner")

    logger.info(f"Variant start position must be within an NMD region.")
    logger.info(f"ClinVar variants with NMD region annotation: {len(df)}")

    # Annotate variants with constraint data
    df = df.merge(constraint, how="left")

    logger.info(f"Variants after merging constraint data: {len(df)}")
    logger.info(f"Non-null values:\n{df.info()}")

    logger.info("Writing to output.")
    df.to_csv(C.CLINVAR_LOF_ANNOTATED, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    df = main()
