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


def gost(query, background, **kwargs):
    kwargs.setdefault("sources", ["GO:MF", "GO:BP", "HP"])
    kwargs.setdefault("no_iea", True)
    kwargs.setdefault("domain_scope", "custom_annotated")
    kwargs.setdefault("significance_threshold_method", "bonferroni")
    kwargs.setdefault("no_evidences", True)

    gp = GProfiler(return_dataframe=True)

    profile = gp.profile(
        query=query,
        background=background,
        **kwargs,
    )  # [["source","native","name","p_value", "query", "highlight"]]

    return profile


def main():
    """Run as script."""

    gene_sets = [pd.read_csv(p).iloc[:, 0].tolist() for p in _PATHS]
    gene_sets = dict(zip(_NAMES, gene_sets))

    background_genes = pd.read_csv(C.GENE_LIST_ALL).iloc[:, 0].tolist()

    profile = gost(gene_sets, background_genes)

    profile["p_value_rank"] = profile.groupby(["query", "source"])["p_value"].rank(
        method="dense"
    )
    profile["enrichment"] = (profile["intersection_size"] / profile["query_size"]) / (
        profile["term_size"] / profile["effective_domain_size"]
    )
    profile["enrichment_rank"] = profile.groupby(["query", "source"])[
        "enrichment"
    ].rank(ascending=False, method="dense")

    profile = profile[
        [
            "query",
            "source",
            "native",
            "name",
            "p_value",
            "p_value_rank",
            "enrichment",
            "enrichment_rank",
        ]
    ]

    # profile = profile[profile["rank"].le(10)]

    return profile


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
