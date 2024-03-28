"""Module docstring."""

# Imports
import logging
from pathlib import Path

import pandas as pd

import src
from src import constants as C
from src.constraint import constraint_statistics

# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"


# Logging
logger = logging.getLogger(__name__)


# Functions
def read_gene_ids(path):
    return pd.read_csv(path, sep="\t", usecols=["gene_id", "transcript_id"]).set_axis(
        ["ensg", "enst"], axis="columns"
    )


def get_constrained_gnomad_genes(gnomad, gene_ids):

    # Find constrained transcripts
    m1 = gnomad.pli > 0.9
    m2 = gnomad.loeuf < 0.6

    gnomad_strong = gnomad[m1 | m2]["enst"].drop_duplicates()

    logger.info(f"gnomAD constrained transcripts: {len(gnomad_strong)}")

    # Get gene ids
    gnomad_strong = gene_ids.merge(gnomad_strong, how="inner", validate="one_to_one")[
        "ensg"
    ]
    assert gnomad_strong.duplicated().sum() == 0, "Duplicated ENSG IDs."

    logger.info(f"gnomAD constrained gene ids: {len(gnomad_strong)}")

    return gnomad_strong


def main():
    """Run as script."""

    # Read data
    gene_ids = read_gene_ids(C.CANONICAL_CDS_GENE_IDS)
    gnomad = constraint_statistics.get_gnomad_constraint(C.GNOMAD_V4_CONSTRAINT)

    # Define gene lists
    ## All genes with a protein-coding canonical transcript
    _all = gene_ids["ensg"].drop_duplicates()
    logger.info(f"All genes: {len(_all)}")

    ## Constrained genes in gnomAD
    gnomad_strong = get_constrained_gnomad_genes(gnomad, gene_ids)



    # Write to output
    _all.to_csv(C.GENE_LIST_ALL, index=False, header=None)
    gnomad_strong.to_csv(C.GENE_LIST_GNOMAD_CST, index=False, header=None)

    return gene_ids, gnomad, gnomad_strong


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
