"""Module docstring."""

import logging
from pathlib import Path

import pandas as pd
from gprofiler import GProfiler

import src
from src import constants as C

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_CONSTRAINED_PATHS = [
    C.GENE_LIST_GNOMAD_CST,
    C.GENE_LIST_NMD_TARGET,
    C.GENE_LIST_START_PROX,
    C.GENE_LIST_LONG_EXON,
    C.GENE_LIST_DISTAL,
]
_CONSTRAINED_NAMES = ["gnomAD", "NMD target", "Start proximal", "Long exon", "Distal"]


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
    )

    return profile


def get_enrichment_stats(profile, keep_rank=10):
    """Quantify and rank relative enrichment of terms."""

    profile["enrichment"] = (profile["intersection_size"] / profile["query_size"]) / (
        profile["term_size"] / profile["effective_domain_size"]
    )

    profile["enrichment_rank"] = profile.groupby(["query", "source"])[
        "enrichment"
    ].rank(ascending=False, method="first")

    # Keep only top-ranking terms
    profile = profile[profile["enrichment_rank"].le(keep_rank)]

    return profile


def tidy_data(profile, background="all_genes"):
    profile = profile.assign(background=background)[
        [
            "query",
            "background",
            "source",
            "native",
            "name",
            "enrichment",
            "enrichment_rank",
            "p_value",
        ]
    ].sort_values(["source", "query", "enrichment_rank"])

    return profile


def gene_set_enrichment():
    return None


def main():
    """Run as script."""

    # Read gene lists
    list_genes = lambda x: pd.read_csv(x).iloc[:, 0].tolist()

    _all = list_genes(C.GENE_LIST_ALL)
    gnomad = list_genes(C.GENE_LIST_GNOMAD_CST)
    target = list_genes(C.GENE_LIST_NMD_TARGET)
    start = list_genes(C.GENE_LIST_START_PROX)
    long_exon = list_genes(C.GENE_LIST_LONG_EXON)
    distal = list_genes(C.GENE_LIST_DISTAL)

    # Gene set enrichment: constrained genes vs all genes
    profile = gost(distal, _all)#.pipe(get_enrichment_stats).pipe(tidy_data)

    # Regional constraint vs gnomAD
    # regional_gene_sets = [pd.read_csv(p).iloc[:, 0].tolist() for p in _REGIONAL_PATHS]
    # regional_gene_sets = dict(zip(_REGIONAL_NAMES, regional_gene_sets))

    # gnomad_constrained_genes = pd.read_csv(C.GENE_LIST_GNOMAD_CST).iloc[:, 0].tolist()

    # region_profile = (
    #     gost(regional_gene_sets, gnomad_constrained_genes)
    #     .pipe(get_enrichment_stats)
    #     .pipe(tidy_data, "gnomad")
    # )

    return profile


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()

# Don't bother with under-representation analysis - the numbers are too low.
