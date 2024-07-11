"""Simplify the parsed genemap2 annotation.

Notes:
- Non-disease phenotypes (e.g. provisional, or susceptibility phenotypes) are dropped.
- OMIM entries where inheritance information is missing are dropped.
- The inheritance patterns from OMIM are simplified here.
- Several genes have duplicate entries in OMIM (multiple phenotypes and / or inheritance
  modes). These genes are also duplicated in the merged annotation.
"""

import logging
from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/interim/genemap2_parsed.tsv"
_FILE_OUT = "data/interim/genemap2_simple.tsv"

logger = logging.getLogger(__name__)


def read_parsed_genemap_data(path):
    """Read parsed genemap data to memory."""

    omim = pd.read_csv(
        path,
        sep="\t",
        usecols=["ensg", "phenotype", "inheritance"],
    )

    logger.info(f"Entries in parsed genemap data: {len(omim)}")

    omim = omim.dropna(subset=["ensg", "inheritance"]).drop_duplicates()

    logger.info(f"Entries after dropping missing values and duplicates: {len(omim)}")

    return omim


def exclude_non_disease_phenotypes(omim):
    """Exclude susceptibility, provisional, and non disease phenotypes."""

    m1 = omim.phenotype.str.startswith("[")  # Non-disease
    m2 = omim.phenotype.str.startswith("{")  # Susceptibility
    m3 = omim.phenotype.str.startswith("?")  # Provisional

    logger.info(f"Non-disease phenotypes:{m1.sum()}")
    logger.info(f"Susceptibility phenotypes: {m2.sum()}")
    logger.info(f"Provisional phenotypes: {m3.sum()}")

    omim = omim[~(m1 | m2 | m3)]

    logger.info(f"Remaining entries: {len(omim)}")

    return omim


def sanitise_inheritance_modes(omim):
    """Categorise as AD, AR, X-linked, or Other."""

    logger.info(f"Inheritance mode value counts:\n{omim.inheritance.value_counts()}")

    m1 = omim.inheritance.str.contains("X-linked")
    m2 = omim.inheritance.str.startswith("Autosomal recessive")
    m3 = omim.inheritance.str.startswith("Autosomal dominant")

    omim.loc[m1, "inheritance"] = "X-linked"
    omim.loc[~(m1 | m2 | m3), "inheritance"] = "Other"

    logger.info(
        f"Duplicated entries after cleaning (these are dropped): {omim.duplicated().sum()}"
    )

    omim = omim.drop_duplicates()

    logger.info(
        f"Inheritance mode value counts after cleaning:\n{omim.inheritance.value_counts()}"
    )
    logger.info(
        f"Unique genes in each inheritance group:\n{omim.groupby('inheritance').ensg.nunique()}"
    )
    logger.info(f"Unique monogenic disease genes: {omim.ensg.nunique()}")

    return omim


def main():
    omim = (
        read_parsed_genemap_data(_FILE_IN)
        .pipe(exclude_non_disease_phenotypes)
        .pipe(sanitise_inheritance_modes)
    )

    logger.info("Writing to output.")
    omim.to_csv(_FILE_OUT, sep="\t", index=False)

    return omim


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()