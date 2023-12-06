"""Merge NMD, phyloP, AlphaMissense, and pext annotations."""

# Imports
from pathlib import Path

import numpy as np
import pandas as pd

from src import constants as C
from src import setup_logger


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
        nrows=100,  #! Testing
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
        nrows=100,  #! Testing
    )

    logger.info(f"phyloP annotations: {len(phylop)}")

    return phylop


def read_pext_annotations(path):
    """Get base level pext annotations."""

    logger.info("Getting pext annotations")

    pext = pd.read_csv(
        path,
        sep="\t",
        nrows=100,  #! Testing
        header=None,
        names=["chr", "start", "pos", "ensg", "pext"],
        usecols=["chr", "pos", "ensg", "pext"],
    )

    logger.info(f"pext annotations: {len(pext)}")

    return pext


def main():
    """Run as script."""
    nmd = read_nmd_annotations(C.NMD_ANNOTATIONS)
    phylop = read_phylop_annotations(C.PHYLOP_CDS_SCORES)
    pext = read_pext_annotations(C.PEXT_BED_38)

    return pext  #! Testing


if __name__ == "__main__":
    main()


# # Read pext data into memory
# pext = (
#     pd.read_csv(
#         "../outputs/pext_38.bed",
#         sep="\t", header=None, names=["chr","start","end","ensg","pext"],
#         usecols=["chr","end","ensg","pext"]
#     )
#     .rename(columns={"end":"pos"})
#     .drop_duplicates()
#     .drop_duplicates(["chr","pos","ensg"], keep=False)
# )
# print(f"Valid pext annotations: {len(pext)}")

# # Read gene and transcript ids into memory
# ids = (
#     pd.read_csv(
#         "../outputs/gene_ids.tsv",
#         sep="\t",
#         header=0,
#         names=["ensg","enst","hgnc"],
#         usecols=["ensg","enst"]
#     )
# )
# ids["ensg"] = ids["ensg"].str.split(".").str[0]
# ids["enst"] = ids["enst"].str.split(".").str[0]

# ids = ids.drop_duplicates()

# pext = pext.merge(ids, how="inner").drop("ensg", axis=1)
# print(f"Valid pext annotations in genes with a MANE transcript: {len(pext)}")


# # ### HMC annotations
# hmc = (
#     pd.read_csv(
#         "../outputs/hmc_38.tsv",
#         sep="\t",
#         header=None,
#         names=["chr","pos","hmc"]
#     )
#     .sort_values(["chr","pos","hmc"])
#     .drop_duplicates(["chr","pos"]) # Keep the lowest HMC score (most constrained) per site
# )
# print(f"Number of HMC annotations: {len(hmc)}")


# # ### Constraint annotations

# # Read the constraint data into memory
# constraint = (
#     pd.read_csv(
#         "../outputs/expected_variants_all_regions_no_cpg_stats.tsv",
#         sep="\t",
#         usecols=["region", "enst", "csq", "n_obs", "oe", "z", "p", "fdr_p"],
#     )
#     .pivot( # We need, for example, synonymous Z-scores for later filtering
#         index=["region", "enst"],
#         columns="csq",
#         values=["n_obs", "oe", "z", "p", "fdr_p"],
#     )
#     .swaplevel(
#         axis=1,
#     )
#     .reset_index(
#         drop=False,
#     )
# )


# # Find constrained and unconstrained regions

# ## The columns are a multi-index which need to be merged
# constraint.columns = ["_".join(x).strip("_") for x in constraint.columns.values]

# ## Keep only the relevant columns
# constraint = constraint[
#     [
#         "region",
#         "enst",
#         "nonsense_n_obs",
#         "nonsense_oe",
#         "synonymous_z",
#         "nonsense_p",
#         "nonsense_fdr_p",
#     ]
# ]

# ## Filter for constrained and unconstrained regions / transcripts
# m1 = constraint["nonsense_oe"] < 0.35
# m2 = constraint["synonymous_z"] > -1
# m3 = constraint["nonsense_fdr_p"] < 0.05

# m4 = constraint["nonsense_p"] >= 0.05
# m5 = constraint["nonsense_n_obs"] >= 1

# constraint.loc[m1 & m2 & m3, "constraint"] = "constrained"
# constraint.loc[m4 & m5, "constraint"] = "unconstrained"

# ## Drop irrelevant columns
# constraint = constraint[["region", "enst", "constraint"]]

# ## Print the counts of constrained and unconstrained regions
# print(constraint.groupby(["region"])["constraint"].value_counts())


# # ## Merge annotations


# # NMD and phyloP
# df = nmd.merge(phylop, how="left")
# print(f"Sites after merging NMD and phyloP annotations: {len(df)}")
# print(f"Sites with a phyloP annotation: {len(df) - df.phylop.isna().sum()}")

# # pext
# df = df.merge(pext, how="left")
# print(f"Sites after merging pext annotations: {len(df)}")
# print(f"Sites with a pext annotation: {len(df) - df.pext.isna().sum()}")

# # hmc
# df = df.merge(hmc, how="left")
# print(f"Sites after merging with HMC annotation: {len(df)}")
# print(f"Sites with an HMC annotation: {len(df) - df.hmc.isna().sum()}")


# # In order to get transcript-level statistics, we copy the dataframe and overwrite the "region" annotation.
# _ = df.copy().assign(region="transcript")
# df = pd.concat([df, _])


# # Merge with constraint annotations
# df = df.merge(constraint, how="inner")


# # ## Statistics


# def se(p, n):
#     """Calculate the standard error of a propotion."""
#     return np.sqrt((p * (1 - p))/n)


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
