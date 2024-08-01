"""Annotate sites with CADD scores."""

# Imports
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import src
from src.orthogonal_metrics import merge_orthogonal_annotations as moa
from src.data import observed_variants

_CADD = "data/interim/cadd_scores_coding.tsv"
_NMD = "data/interim/nmd_annotations.tsv"
_VEP = "data/interim/cds_all_possible_snvs_vep_tidy.tsv"
_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
_FILE_OUT = "data/interim/cadd_scores_coding_annotated.tsv"
_CADD_HEADER = "chr pos ref alt cadd_raw cadd_phred".split()
_CADD_USECOLS = "chr pos ref alt cadd_phred".split()
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_DTYPES = {
    "pos": np.int32,
    "cadd_phred": np.float16,
}


# Logging
logger = logging.getLogger(__name__)


# Functions
def read_cadd(path, **kwargs):
    """Read CADD data."""

    logger.info("Reading CADD data.")

    cadd = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=_CADD_HEADER,
        usecols=_CADD_USECOLS,
        low_memory=False,
        dtype=_DTYPES,
        **kwargs,
    )

    logger.info(f"CADD annotations: {len(cadd)}")

    return cadd


def tidy_cadd(df):
    """Tidy the CADD annotations for later merging."""

    logger.info("Dropping CADD duplicates.")
    df = df.drop_duplicates().copy()
    logger.info(f"CADD annotations (duplicates dropped): {len(df)}")

    logger.info('Adding "chr" prefix to CADD contig names.')
    df["chr"] = ["".join(["chr", str(x)]) for x in df["chr"]]

    return df


def main():
    """Run as script."""

    # Read the data
    cadd = read_cadd(_CADD).pipe(tidy_cadd)
    vep = observed_variants.get_vep_annotations(_VEP)
    nmd = moa.read_nmd_annotations(_NMD, dtype=_DTYPES)
    constraint = moa.read_regional_nonsense_constraint(_CONSTRAINT)

    # Merge annotations
    logger.info("Merging VEP and CADD annotations.")
    df = vep.merge(cadd, how="inner", validate="many_to_one")
    del vep, cadd  # Manage memory
    logger.info(f"Remaining variants: {len(df)}")

    logger.info("Merging NMD annotations.")
    df = df.merge(nmd, how="inner", validate="many_to_one")
    del nmd  # Manage memory
    logger.info(f"Remaining variants: {len(df)}")

    logger.info("Merging constraint annotations.")
    df = df.merge(constraint, how="inner", validate="many_to_one")
    logger.info(f"Remaining variants: {len(df)}")

    # Write to output
    logger.info("Writing to output.")
    df.to_csv(_FILE_OUT, sep="\t", index=False)

    logger.info("Done.")

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
