# # Observed variants
# This script identifies the observed and possible variants in UKB.
#
# Variants are grouped by variant context, transcript, or NMD-region. The script aggregates the number observed, the number possible, and the mean mutability for each grouping. The summary data are saved to a .tsv outputs.

# Imports
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as _stats
import seaborn as sns
import statsmodels.formula.api as smf
from statsmodels.stats.proportion import proportions_ztest

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
def get_observed_variants(path):
    """Get observed variants in gnomAD v4.0"""

    obs = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=_VCF_HEADER + ["ac", "an"],
        usecols=["chr", "pos", "ref", "alt", "ac", "an"],
        dtype=_DATATYPES,
    ).assign(obs=1)

    return obs


def get_vep_annotations(path):
    """Retreive VEP annotations of all possible SNVs."""

    vep = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=_VCF_HEADER,
        dtype=_DATATYPES,
        usecols=["chr", "pos", "ref", "alt", "info"],
    )

    return vep

def get_methylation_data(path):
    """Get gnomAD non-coding methylation data."""
    meth = pd.read_csv(
        path,
        "../data/grch38_cpg_methylation.tsv",
        sep="\t",
        header=0,
        names=["ix", "chr", "pos", "alleles", "lvl"],
        usecols=["chr", "pos", "lvl"],
    )

    return meth


def main():
    """Run the script."""
    tri = pd.read_csv(C.CDS_ALL_SNVS_TRI_CONTEXT, sep="\t", dtype=_DATATYPES)

    return df # TODO Testing

# gnomAD methylation data



# Mutation rates
mu = pd.read_csv(
    "../data/gnomad_nc_mutation_rates.tsv",
    sep="\t",
    names=[
        "tri",
        "ref",
        "alt",
        "lvl",
        "variant_type",
        "mu",
        "pos",
        "obs",
        "po",
        "ppo",
    ],
    header=0,
    usecols=["tri", "ref", "alt", "lvl", "mu", "variant_type"],
)

# Simplify variant type annotations
mu = mu.replace(
    {
        "transversion": "non-CpG",
        "non-CpG transition": "non-CpG",
    }
)

# Mutation rates are only available for 32 codons. We need to reverse-complement for the remainder.
complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

# Replace ref and alt alleles
_mu = mu.copy().replace(complement)

# Reverse-complement trinucleotide contexts
_mu["tri"] = pd.Series(["".join([complement[y] for y in x])[::-1] for x in mu.tri])

# Merge original and reverse-complemented data
mu = pd.concat([mu, _mu])

# NMD data
nmd = pd.read_csv(
    "../data/nmd_annotations.tsv",
    sep="\t",
    usecols=["chr", "pos", "transcript_id", "nmd_definitive"],
).rename(columns={"transcript_id": "enst", "nmd_definitive": "nmd"})

# ## Merge annotations

# Merge VEP, context, NMD, and observed variant annotations
df = vep.merge(tri, how="left")
df = df.merge(nmd, how="left")
df = df.merge(obs, how="left").fillna(0)

# Merge methylation annotations
variant_types = mu[["tri", "ref", "alt", "variant_type"]].drop_duplicates()

df = df.merge(variant_types, how="left")
df = df.merge(meth, how="left")

# All non-CpG sites have lvl 0
df.loc[df["variant_type"] != "CpG", "lvl"] = 0
df.lvl = df.lvl.astype(int)
df.obs = df.obs.astype(int)

# Merge with mutability data
df = df.merge(mu, how="left")


# ## Summarise the data


# ### Rare synonymous variants for expectation model


# Subset to synonymous variants only
syn = df[df["csq"] == "synonymous"]

# Mask contexts in which a synonymous variant is generally not possible.
# (NB synonymous variants in these contexts can only occur at exon-intron junctions)
m1 = (syn.tri == "AGT") & ((syn.alt == "C") | (syn.alt == "T"))
m2 = (syn.tri == "AAT") & ((syn.alt == "C") | (syn.alt == "T"))
m3 = (syn.tri == "ACT") & ((syn.alt == "G") | (syn.alt == "A"))
m4 = (syn.tri == "ATT") & ((syn.alt == "G") | (syn.alt == "A"))

# Mask rare variants
m5 = syn["ac"] == 0
m6 = (syn["ac"] / syn["an"]) < 0.001

# Apply filters
syn = syn[~(m1 | m2 | m3 | m4) & (m5 | m6)]

# Group by variant context
dfg = (
    syn.groupby(["tri", "ref", "alt", "variant_type", "lvl"])
    .agg({"mu": "mean", "obs": "sum", "pos": "count"})
    .reset_index()
)

# Write to output
dfg.to_csv("../outputs/observed_variants_stats_synonymous.tsv", sep="\t", index=False)


# ### Observed variants by region


# Transcripts
dfg = (
    df.groupby(["enst", "csq", "variant_type"])
    .agg(n_pos=("pos", "count"), n_obs=("obs", "sum"), mu=("mu", "mean"))
    .reset_index()
)

# Save to output
dfg.to_csv("../outputs/observed_variants_stats_transcript.tsv", sep="\t", index=False)


# NMD regions
dfg = (
    df.groupby(["enst", "csq", "variant_type", "nmd"])
    .agg(n_pos=("pos", "count"), n_obs=("obs", "sum"), mu=("mu", "mean"))
    .reset_index()
)

# Save to output
dfg.to_csv("../outputs/observed_variants_stats_nmd.tsv", sep="\t", index=False)


# # Upload outputs to UKB RAP
# dx upload --destination outputs/ ../outputs/observed_variants_stats_transcript.tsv ../outputs/observed_variants_stats_nmd.tsv
