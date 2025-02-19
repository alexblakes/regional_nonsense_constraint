"""Merge NMD, phyloP, AlphaMissense, and pext annotations."""



import pandas as pd

import src


_NMD = "data/interim/nmd_annotations.tsv"
_PHYLOP = "data/interim/phylop_cds_sites.tsv"
_PEXT = "data/interim/pext_38.bed"
_GENE_IDS = "data/interim/gene_ids.tsv"
_ALPHA_MIS = "data/interim/alpha_missense_tidy.tsv"
_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
_FILE_OUT = "data/interim/cds_sites_phylop_pext_missense.tsv"
_STATS_OUT = "data/final/phylop_pext_missense_annotations_stats.tsv"
_METRICS = ["phylop", "pext", "alpha_mis"]

logger = src.logger


def read_nmd_annotations(path, **kwargs):
    """Get per-site NMD annotations."""

    logger.info("Getting NMD annotations.")

    nmd = pd.read_csv(
        path,
        sep="\t",
        usecols=["chr", "pos", "transcript_id", "nmd_definitive"],
        **kwargs,
    ).rename(columns={"nmd_definitive": "region", "transcript_id": "enst"})

    logger.info(f"NMD annotations: {len(nmd)}")
    logger.info(
        f"Duplicated by site / enst: {nmd.duplicated(['chr','pos','enst']).sum()}"
    )

    return nmd


def read_phylop_annotations(path):
    """Get per-site phyloP annotations."""

    logger.info("Getting phyloP annotations.")

    phylop = pd.read_csv(
        path,
        sep="\t",
        # nrows=10000,  #! Testing
    )

    logger.info(f"phyloP annotations: {len(phylop)}")

    return phylop


def read_pext_annotations(path):
    """Get base level pext annotations."""

    logger.info("Getting pext annotations")

    pext = pd.read_csv(
        path,
        sep="\t",
        # nrows=10000,  #! Testing
        header=None,
        names=["chr", "start", "pos", "ensg", "pext"],
        usecols=["chr", "pos", "ensg", "pext"],
    )

    logger.info(f"pext annotations: {len(pext)}")

    return pext


def read_gene_ids(path):
    """Get gene / transcript IDs."""

    ids = pd.read_csv(
        path,
        sep="\t",
        usecols=["gene_id", "transcript_id"],
    ).set_axis(["ensg", "enst"], axis=1)

    return ids


def read_alpha_missense(path):
    """Get tidied AlphaMissense scores."""

    logger.info("Getting AlphaMissense annotations.")

    am = pd.read_csv(
        path,
        sep="\t",
        # nrows=10000,  #! Testing
    ).rename(columns={"alpha_missense_min": "alpha_mis"})

    logger.info(f"AlphaMissense annotations: {len(am)}")

    return am


def read_regional_nonsense_constraint(path, **kwargs):
    """Get regional nonsense constraint annotations."""

    constraint = pd.read_csv(
        path,
        sep="\t",
        usecols=["enst", "region", "constraint"],
        **kwargs,
    )

    return constraint


def get_mean_score(df, metric, group_by_constraint=False):
    """Find the mean score of the annotation in NMD regions."""

    if not group_by_constraint:
        df = df.groupby("region")
    elif group_by_constraint:
        df = df.groupby(["region", "constraint"], dropna=False)

    df = (
        df[metric]
        .agg(["mean", "count", "sem"])
        .set_axis([metric, "n", "sem"], axis=1)
        .assign(ci_95=lambda x: x["sem"] * 1.96)
        .reset_index(drop=False)
    )

    if not group_by_constraint:
        df = df.assign(constraint="all")

    df = df.melt(
        id_vars=["region", "n", "sem", "ci_95", "constraint"],
        var_name="metric",
        value_name="mean",
    )[["region", "constraint", "metric", "mean", "n", "sem", "ci_95"]].sort_values(
        ["constraint", "region"]
    )

    return df


def combine_mean_scores_across_groups(df, annotation):
    """Concatenate mean scores for regions and constraint groups."""
    return pd.concat(
        [
            get_mean_score(df, annotation, False),
            get_mean_score(df, annotation, True),
        ]
    )


def main():
    """Run as script."""

    nmd = read_nmd_annotations(_NMD)
    phylop = read_phylop_annotations(_PHYLOP)
    pext = read_pext_annotations(_PEXT)
    ids = read_gene_ids(_GENE_IDS)
    am = read_alpha_missense(_ALPHA_MIS)
    constraint = read_regional_nonsense_constraint(_CONSTRAINT)

    # Merge the annotations

    ## PhyloP
    logger.info("Merging NMD and phylop annotations.")
    df = nmd.merge(phylop, how="left")
    logger.info(f"Merged NMD and phyloP annotations: {len(df)}")
    logger.info(f"Sites with a phyloP annotation: {len(df) - df.phylop.isna().sum()}")

    ## pext
    logger.info("Getting transcript IDs for pext annotations.")
    pext = pext.merge(ids, how="inner").drop("ensg", axis=1)
    logger.info(f"Remaining pext annotations after getting transcript IDs: {len(pext)}")

    logger.info("Merging pext annotations.")
    df = df.merge(pext, how="left")
    logger.info(f"Sites after merging pext annotations: {len(df)}")
    logger.info(f"Sites with a pext annotation: {len(df) - df.pext.isna().sum()}")

    ## AlphaMissense
    logger.info("Merging AlphaMissense annotations.")
    df = df.merge(am, how="left")
    logger.info(f"Sites after merging with AlphaMissense scores: {len(df)}")
    logger.info(
        f"Sites with an AlphaMissense score: {len(df) - df.alpha_mis.isna().sum()}"
    )

    ## Constraint annotations
    ### Get transcript-level annotations
    logger.info("Adding 'transcript' to the set of regions.")
    df = pd.concat([df, df.copy().assign(region="transcript")])
    logger.info(f"Sites after adding 'transcript' to the set of regions: {len(df)}")

    logger.info("Merging constraint annotations.")
    df = df.merge(constraint, how="inner")
    logger.info(f"Sites after merging with constraint annotations: {len(df)}")
    logger.info(
        f"Sites per region and constraint annotation:\n{df.groupby('region').constraint.value_counts(dropna=False)}"
    )

    # Mean scores for each annotation
    stats = pd.concat([combine_mean_scores_across_groups(df, x) for x in _METRICS])

    # Write to output
    logger.info("Writing to output.")
    df.to_csv(_FILE_OUT, sep="\t", index=False)
    stats.to_csv(_STATS_OUT, sep="\t", index=False)

    return stats  #! Testing


if __name__ == "__main__":
    src.add_log_handlers()
    main()
