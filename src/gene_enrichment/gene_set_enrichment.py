"""Module docstring."""

import logging
from pathlib import Path

import pandas as pd
from gprofiler import GProfiler

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_PATHS = [
    C.GENE_LIST_GNOMAD_CST,
    C.GENE_LIST_NMD_TARGET,
    C.GENE_LIST_START_PROX,
    C.GENE_LIST_LONG_EXON,
    C.GENE_LIST_DISTAL,
]
_NAMES = ["gnomAD", "NMD target", "Start proximal", "Long exon", "Distal"]

logger = logging.getLogger(__name__)


def main():
    """Run as script."""

    gene_sets = [pd.read_csv(p).iloc[:,0].tolist() for p in _PATHS]
    gene_sets = dict(zip(_NAMES, gene_sets))

    background_genes = pd.read_csv(C.GENE_LIST_ALL).iloc[:,0].tolist()

    gp = GProfiler(return_dataframe=True)

    profile = gp.profile(
        organism="hsapiens",
        query=gene_sets,
        sources=["GO:MF", "GO:BP", "HP"],
        user_threshold=0.05,
        combined=True,
        no_iea=True,
        domain_scope="custom_annotated",
        significance_threshold_method="g_SCS",
        background=background_genes,
        no_evidences=True,
        # highlight=True,
    )#[["source","native","name","p_values"]]

    return profile


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
