"""
Annotate observed and possible variants.

Merge annotations for trinucleotide context, vep consequences, NMD annotations, 
methylation, and mutability for all observed and possible SNVs.
"""

# Imports
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from src import setup_logger
from src import constants as C

# Module constants
_VCF_HEADER = ["chr", "pos", "id", "ref", "alt", "qual", "filter", "info"]
_DATATYPES = defaultdict(lambda: "str").update(
    {"pos": np.int32, "ac": np.int32, "an": np.int32}
)


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_trinucleotide_contexts(path):
    """Get trinucleotide contexts."""

    logger.info("Getting trinucleotide sequence contexts.")

    tri = pd.read_csv(
        path,
        sep="\t",
        dtype=_DATATYPES,
        # nrows=10000,  #! Testing
    )

    logger.info(f"Positions with trinucleotide sequence contexts: {len(tri)}")

    return tri


def get_vep_annotations(path):
    """Retreive VEP annotations of all possible SNVs."""

    logger.info("Getting VEP annotations for all possible SNVs.")

    vep = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=["chr", "pos", "ref", "alt", "csq", "enst"],
        dtype=_DATATYPES,
        # nrows=10000  # ! Testing
    )

    logger.info(f"VEP-annotated possible SNVs: {len(vep)}")

    return vep


def get_nmd_annotations(path):
    """Get NMD annotations for all CDS sites."""

    logger.info("Getting NMD annotations for all CDS sites.")

    nmd = pd.read_csv(
        path,
        sep="\t",
        usecols=["chr", "pos", "transcript_id", "nmd_definitive"],
        dtype=_DATATYPES,
        # nrows=10000,  # ! Testing
    ).rename(columns={"transcript_id": "enst", "nmd_definitive": "nmd"})

    logger.info(f"NMD annotations: {len(nmd)}")

    return nmd


def get_observed_variants(path):
    """Get observed variants in gnomAD v4.0"""

    logger.info("Getting observed variants in gnomAD v4.0.")

    obs = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=_VCF_HEADER + ["ac", "an"],
        usecols=["chr", "pos", "ref", "alt", "ac", "an"],
        dtype=_DATATYPES,
        # nrows=10000,  #! Testing
    ).assign(obs=True)

    logger.info(f"Observed variants: {len(obs)}")
    logger.info(
        f"Observed variant duplicates: {obs.duplicated(['chr','pos', 'ref','alt']).sum()}"
    )

    return obs


def get_methylation_data(path):
    """Get methylation annotations for all CpG sites."""

    logger.info("Getting methylation annotations for all CpG sites.")

    meth = pd.read_csv(
        path,
        sep="\t",
        header=0,
        names=["ix", "chr", "pos", "alleles", "lvl"],
        usecols=["chr", "pos", "lvl"],
        # nrows=10000,  # ! Testing
    )

    logger.info(f"CpG methylation annotations: {len(meth)}")

    return meth


def assign_methylation_levels(df):
    """Assign methylation levels to non-CpG variants."""

    logger.info(f"Missing lvl annotations: {df['lvl'].isna().sum()}")

    # All non-CpG sites have methylation level 0.
    df.loc[df["variant_type"] != "CpG", "lvl"] = 0
    df.lvl = df.lvl.astype(int)  # For later merging with mutability data

    logger.info(
        f"Missing lvl annotations after assigning non-CpGs to lvl = 0: {df['lvl'].isna().sum()}"
    )

    return df


def get_coverage_data(path):
    """Get median coverage at all sites."""

    logger.info("Getting coverage data.")

    cov = pd.read_csv(
        path,
        sep="\t",
        usecols=["locus", "median_approx"],
        dtype={"locus": "str", "median_approx": "int16"},
        # nrows=10000,  #! Testing
    ).rename(columns={"median_approx": "median_coverage"})

    locus = cov.locus.str.split(":")
    cov["chr"] = locus.str[0]
    cov["pos"] = locus.str[1].astype("int32")
    cov = cov.drop("locus", axis=1)[["chr", "pos", "median_coverage"]].drop_duplicates()

    logger.info(f"Coverage annotations: {len(cov)}")

    return cov


def address_missing_values(df):
    """Address missing values in the data."""

    logger.info(
        f"Variants lacking an nmd annotation are dropped: N = {df.nmd.isna().sum()}"
    )
    logger.info(
        f"Variants lacking a variant_type annotation are dropped:\n{df[df['variant_type'].isna()]}"
    )

    df = df.dropna(subset=["nmd", "variant_type"])

    logger.info(f"Variants after addressing missing values: {len(df)}")
    logger.info(f"Missing values:\n{df.isna().sum()}")

    return df


def drop_chrm_sites(df):
    """Drop chrM sites."""

    m = df[df["chr"] == "chrM"]

    logger.info(f"Unique chromosomes: {df.chr.nunique()}")
    logger.info(f"chrM variants: {len(m)}")
    logger.info(f"chrM transcripts: {m.enst.nunique()}")
    logger.info(f"Observed chrM variants: {m.obs.sum()}")
    logger.info(f"chrM sites are dropped.")

    df = df[df["chr"] != "chrM"]

    logger.info(f"Variants after dropping chrM sites: {len(df)}")

    return df


def log_summary_data(df):
    """Log summary data."""

    # Summary statistics (all variants)
    logger.info(
        f"Duplicated variants: {df.duplicated(['chr','pos','ref','alt', 'enst']).sum()}"
    )
    logger.info(f"Unique transcripts: {df.enst.nunique()}")
    logger.info(f"Consequence value counts:\n{df.csq.value_counts()}")
    logger.info(f"CpG value counts:\n{df.variant_type.value_counts()}")
    logger.info(f"NMD regions value counts:\n{df.nmd.value_counts()}")

    # Summary statistics (CpGs)
    cpg = df[df["variant_type"] == "CpG"]

    logger.info(f"CpG levels:\n{cpg.lvl.value_counts()}")

    # Summary statistics (Observed variants)
    obs = df[df["obs"] == True]

    logger.info(
        f"Observed variants consequence value counts:\n{obs.csq.value_counts()}"
    )


def main():
    """Run the script."""

    # Get datasets
    vep = get_vep_annotations(C.VEP_ALL_SNVS_TIDY)  # chr pos ref alt
    tri = get_trinucleotide_contexts(C.CDS_ALL_SNVS_TRI_CONTEXT)  # chr pos ref alt
    nmd = get_nmd_annotations(C.NMD_ANNOTATIONS)  # chr pos
    obs = get_observed_variants(C.GNOMAD_PASS_SNVS)
    meth = get_methylation_data(C.GNOMAD_NC_METHYLATION)
    mu = pd.read_csv(C.GNOMAD_NC_MUTABILITY_TIDY, sep="\t")
    variant_types = mu[["tri", "ref", "alt", "variant_type"]].drop_duplicates()
    cov = get_coverage_data(C.GNOMAD_COVERAGE)

    # Merge the annotations
    logger.info("Merging annotations.")

    df = (
        vep.merge(tri, how="left")
        .merge(nmd, how="left")
        .merge(obs, how="left")
        .fillna({"obs": False})
        .merge(variant_types, how="left")
        .merge(meth, how="left")
        .pipe(assign_methylation_levels)
        .merge(mu, how="left")
        .merge(cov, how="left")
    )

    logger.info(f"Variants in raw merged data: {len(df)}")

    # Tidy the data
    df = address_missing_values(df).pipe(drop_chrm_sites)

    # Write logs
    log_summary_data(df)

    # Write to output
    df.to_csv(C.ALL_VARIANTS_MERGED_ANNOTATIONS, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    main()
