"""Merge NMD, phyloP, AlphaMissense, and pext annotations."""

# Imports
from locale import D_T_FMT
from pathlib import Path

import numpy as np
import pandas as pd

from src import constants as C
from src import setup_logger
from src.functional_clinical.alpha_missense_tidy import read_alpha_missense


# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_nmd_annotations(path):
    """Get per-site NMD annotations."""

    logger.info("Getting NMD annotations.")

    nmd = pd.read_csv(
        path,
        sep="\t",
        nrows=10000,  #! Testing
        usecols=["chr", "pos", "transcript_id", "nmd_definitive"],
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
        nrows=10000,  #! Testing
    )

    logger.info(f"phyloP annotations: {len(phylop)}")

    return phylop


def read_pext_annotations(path):
    """Get base level pext annotations."""

    logger.info("Getting pext annotations")

    pext = pd.read_csv(
        path,
        sep="\t",
        nrows=10000,  #! Testing
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
    am = pd.read_csv(
        path,
        sep="\t",
        nrows=10000,  #! Testing
    )

    logger.info(f"AlphaMissense annotations: {len(am)}")

    return am


def read_regional_nonsense_constraint(path):
    """Get regional nonsense constraint annotations."""

    constraint = pd.read_csv(path, sep="\t", usecols=["enst", "region", "constraint"])

    return constraint


def se(p, n):
    """Calculate the standard error of a propotion."""
    return np.sqrt((p * (1 - p)) / n)

def find_sites_meeting_condition(df, annotation, condition):
    """Find sites at which the values for an annotation meet the given condition."""

    df = df[["region", "constraint", annotation]].copy().dropna(subset=annotation)
    df[annotation] = condition(df[annotation])

    return df

def get_proportion_meeting_condition(df, group_on_constraint=False):
    """Find the proportion of sites which fulfil the condition."""

    annotation = df.columns.drop(["region","constraint"])[0]

    if not group_on_constraint:
        df = df.groupby("region")
    elif group_on_constraint:
        df = df.groupby(["region", "constraint"])

    df = df[annotation].agg(["mean","count"]).set_axis([annotation,"n"], axis=1).reset_index(drop=False)

    if not group_on_constraint:
        df.assign(constraint="full_transcript")

    return df

def get_95_ci(df):
    

# # Regions, ignoring constraint annotation
# r = (
#     p.groupby("region")
#     .agg(phylop=("phylop","mean"), n=("phylop","count"))
#     .assign(constraint="all")
#     .reset_index(drop=False)
# )

# # Regions and constraint annotation
# rc = (
#     p.groupby(["region","constraint"])
#     .agg(phylop=("phylop","mean"), n=("phylop","count"))
#     .reset_index(drop=False)
# )

# # Combine
# _p = pd.concat([rc,r])
# _p["ci95"] = 1.96 * se(_p["phylop"], _p["n"])

def main():
    """Run as script."""

    nmd = read_nmd_annotations(C.NMD_ANNOTATIONS)
    phylop = read_phylop_annotations(C.PHYLOP_CDS_SCORES)
    pext = read_pext_annotations(C.PEXT_BED_38)
    ids = read_gene_ids(C.CANONICAL_CDS_GENE_IDS)
    am = read_alpha_missense(C.ALPHA_MISSENSE_TIDY)
    constraint = read_regional_nonsense_constraint(C.REGIONAL_NONSENSE_CONSTRAINT)

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
        f"Sites with an AlphaMissense score: {len(df) - df.alpha_missense_min.isna().sum()}"
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

    # Statistics for each annotation
    _p = find_sites_meeting_condition(df, "phylop", lambda x: x > 2.27)
    _p = pd.concat([get_proportion_meeting_condition(_p), get_proportion_meeting_condition(_p, True)])



    return _p  #! Testing


if __name__ == "__main__":
    main()




# # ## Statistics


# # ### phyloP
# # Calculate the proportion of sites in each grouping with phyloP > 2.27


# # Find sites where phyloP > 2.27
# p = df[["region", "phylop", "constraint"]].copy().dropna(subset="phylop")
# p["phylop"] = np.where(p["phylop"] > 2.27, 1, 0)


# # Regions, ignoring constraint annotation
# r = (
#     p.groupby("region")
#     .agg(phylop=("phylop","mean"), n=("phylop","count"))
#     .assign(constraint="all")
#     .reset_index(drop=False)
# )

# # Regions and constraint annotation
# rc = (
#     p.groupby(["region","constraint"])
#     .agg(phylop=("phylop","mean"), n=("phylop","count"))
#     .reset_index(drop=False)
# )

# # Combine
# _p = pd.concat([rc,r])
# _p["ci95"] = 1.96 * se(_p["phylop"], _p["n"])

# # Write out
# _p.to_csv("../outputs/stats_phylop.tsv", sep="\t", index=False)


# # ### HMC
# # Calculate the proportion of sites in each group where HMC < 1


# # Find sites where HMC < 1
# h = df[["region", "hmc", "constraint"]].copy().dropna(subset="hmc")
# h["hmc"] = np.where(h["hmc"] < 1, 1, 0)


# # Regions, ignoring constraint annotation
# r = (
#     h.groupby("region")
#     .agg(hmc=("hmc","mean"), n=("hmc","count"))
#     .assign(constraint="all")
#     .reset_index(drop=False)
# )

# # Regions and constraint annotation
# rc = (
#     h.groupby(["region","constraint"])
#     .agg(hmc=("hmc","mean"), n=("hmc","count"))
#     .reset_index(drop=False)
# )

# # Combine
# _h = pd.concat([rc,r])
# _h["ci95"] = 1.96 * se(_h["hmc"], _h["n"])

# # Write out
# _h.to_csv("../outputs/stats_hmc.tsv", sep="\t", index=False)


# # ### pext
# # Calculate the mean pext of sites in each group


# # Find sites with a pext annotation
# x = df[["region", "pext", "constraint"]].copy().dropna(subset="pext")


# # Regions, ignoring constraint annotation
# r = (
#     x.groupby("region")
#     .agg(pext=("pext", "mean"), n=("pext","count"), sem=("pext", "sem"))
#     .assign(constraint="all")
#     .reset_index(drop=False)
# )

# # Regions and constraint annotation
# rc = (
#     x.groupby(["region","constraint"])
#     .agg(pext=("pext", "mean"), n=("pext","count"), sem=("pext", "sem"))
#     .reset_index(drop=False)
# )

# # Combine
# _x = pd.concat([rc,r])
# _x["ci95"] = 1.96 * _x["sem"]
# _x = _x.drop("sem", axis=1)

# # Write out
# _x.to_csv("../outputs/stats_pext.tsv", sep="\t", index=False)
