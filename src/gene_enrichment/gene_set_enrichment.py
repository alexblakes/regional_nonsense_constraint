"""Gene set enrichment analysis."""



import pandas as pd
from gprofiler import GProfiler

import src
from src import constants as C



logger = src.logger


def gost(query, background, query_name, bg_name, **kwargs):
    """Run gene set enrichment analysis with new defaults."""

    kwargs.setdefault("sources", ["GO:MF", "GO:BP", "HP"])
    kwargs.setdefault("no_iea", True)
    kwargs.setdefault("domain_scope", "custom_annotated")
    kwargs.setdefault("significance_threshold_method", "bonferroni")
    kwargs.setdefault("no_evidences", True)

    logger.info(
        f"Running GOST with '{query_name}' as query and '{bg_name}' as background."
    )

    gp = GProfiler(return_dataframe=True)

    profile = gp.profile(
        query=query,
        background=background,
        **kwargs,
    )

    logger.info(f"Enriched terms:\n{profile.groupby('source').size()}")

    # Optionally name the background and query gene lists
    if bg_name:
        profile = profile.assign(background=bg_name)
    if query_name:
        profile = profile.assign(query=query_name)

    return profile


def get_enrichment_stats(profile):
    """Quantify and rank the enrichment of terms."""

    profile["enrichment"] = (profile["intersection_size"] / profile["query_size"]) / (
        profile["term_size"] / profile["effective_domain_size"]
    )

    profile["enrichment_rank"] = profile.groupby(["query", "source"])[
        "enrichment"
    ].rank(ascending=False, method="first")

    return profile


def tidy_data(profile, keep_rank=10):
    # Keep only top-ranking terms
    profile = profile[profile["enrichment_rank"].le(keep_rank)]

    # Drop redundant columns
    profile = profile[
        [
            "background",
            "query",
            "source",
            "native",
            "name",
            "enrichment",
            "enrichment_rank",
            "p_value",
        ]
    ].sort_values(["background", "source", "query", "enrichment_rank"])

    return profile


def gene_set_enrichment(query, bg, q_name=None, bg_name=None, **kwargs):
    """Run the gene set enrichment pipeline."""
    return (
        gost(query, bg, q_name, bg_name)
        .pipe(get_enrichment_stats)
        .pipe(tidy_data, **kwargs)
    )


def main():
    """Run as script."""

    # Read gene lists
    list_genes = lambda x: pd.read_csv(x, header=None).iloc[:, 0].tolist()

    _all = list_genes(C.GENE_LIST_ALL)
    gnomad = list_genes(C.GENE_LIST_GNOMAD_CST)
    target = list_genes(C.GENE_LIST_NMD_TARGET)
    start = list_genes(C.GENE_LIST_START_PROX)
    long_exon = list_genes(C.GENE_LIST_LONG_EXON)
    distal = list_genes(C.GENE_LIST_DISTAL)

    # Test constrained genes vs all genes
    query_lists = [gnomad, target, start, long_exon, distal]
    names = ["gnomAD", "NMD target", "Start proximal", "Long exon", "Distal"]

    vs_all = [
        gene_set_enrichment(q, _all, n, "All genes") for q, n in zip(query_lists, names)
    ]
    vs_all = pd.concat(vs_all)

    # Test regionally constrained genes vs gnomAD
    query_lists = [target, start, long_exon, distal]
    names = ["NMD target", "Start proximal", "Long exon", "Distal"]

    vs_gnomad = []
    for q, n in zip(query_lists, names):
        vs_gnomad.append(gene_set_enrichment(q, gnomad, n, "gnomAD"))
    vs_gnomad = pd.concat(vs_gnomad)

    # Combine results
    results = pd.concat([vs_all, vs_gnomad])

    # Write to output
    results.to_csv(C.STATS_GENE_SET_ENRICHMENT, sep="\t", index=False)

    return results


if __name__ == "__main__":
    src.add_log_handlers()
    main()
