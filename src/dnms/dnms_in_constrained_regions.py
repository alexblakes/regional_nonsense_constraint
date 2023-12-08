""" 
Docstring.
"""

# Imports
from pathlib import Path

import pandas as pd

from src import constants as C
from src import setup_logger
from src.functional_clinical import merge_orthogonal_annotations as moa


# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_dnm_data(path):
    """Read DNM data exported from GEL to memory."""

    dnms = (
        pd.read_csv(path, sep="\t")
        .drop(["n_truncating", "constraint"], axis=1)
        .rename(columns={"nmd": "region"})
    )

    logger.info(f"Number of DNMs: {len(dnms)}")

    return dnms


def condense_phenotypes_and_inheritance(omim):
    """Condense OMIM annotations, so each gene is represented only once."""

    omim = (
        omim.groupby("ensg")[["phenotype", "inheritance"]]
        .agg(lambda x: " | ".join(x))
        .drop_duplicates()
        .reset_index(drop=False)
    )

    logger.info(
        f"OMIM entries after condensing phenotype and inheritance annotations: {len(omim)}"
    )
    logger.info(f"Unique ensg IDs: {omim.ensg.nunique()}")

    return omim


def main():
    """Run the script."""

    dnms = read_dnm_data(C.DNMS)
    omim = pd.read_csv(C.OMIM_GENEMAP_SIMPLE, sep="\t").pipe(
        condense_phenotypes_and_inheritance
    )
    constraint = moa.read_regional_nonsense_constraint(C.REGIONAL_NONSENSE_CONSTRAINT)

    # Merge datasets
    df = dnms.merge(omim, how="left")
    logger.info(f"DNMs after merging with OMIM annotation: {len(df)}")

    df = df.merge(constraint, how="left")
    logger.info(f"DNMs after merging with constraint annotation: {len(df)}")

    logger.info(f"Missing values:\n{df.isna().sum()}")

    # Write to output
    logger.info("Writing to output.")
    df.to_csv(C.DNMS_ANNOTATED, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    main()
