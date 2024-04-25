"""Add constraint annotation to ClinVar variants."""

import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_CLINVAR = "data/interim/clinvar_variants_vep_tidy.tsv"
_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
_FILE_OUT = "data/interim/clinvar_variants_constraint.tsv"

logger = logging.getLogger(__name__)


def read_clinvar(path):
    clinvar = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "pos", "enst", "ref", "alt", "csq", "region", "acmg"],
    ).drop_duplicates()

    logger.info(f"ClinVar variants: {len(clinvar)}")
    logger.info(
        f"Duplicated on chr/pos/ref/alt/enst: "
        f"{clinvar.duplicated(['chr','pos','ref','alt','enst']).sum()}"
    )
    logger.info(f"NaNs:\n{clinvar.isna().sum()}")

    return clinvar


def read_constraint(path):
    constraint = pd.read_csv(
        path,
        sep="\t",
        usecols=["enst", "region", "constraint"],
    )

    logger.info(f"Regions in constraint annotation: {len(constraint)}")
    logger.info(
        f"Duplicated on enst/region: "
        f'{constraint.duplicated(["enst","region"]).sum()}'
    )
    logger.info(
        f"Constraint annotation value counts:\n"
        f"{constraint.constraint.value_counts(dropna=False)}"
    )

    return constraint


def main():
    """Run as script."""

    clinvar = read_clinvar(_CLINVAR)
    constraint = read_constraint(_CONSTRAINT)

    df = clinvar.merge(constraint, how="inner")
    logger.info(f"Remaining variants after merge: {len(df)}")

    df.to_csv(_FILE_OUT, sep="\t", index=False)

    return df


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
