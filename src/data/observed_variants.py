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


# Logging
logger = setup_logger(Path(__file__).stem)

# Functions

# ### Download datasets from UKB RAP

# %%
%%bash
dx download \
    -f \
    -o ../data/ \
    data/cds_trinucleotide_contexts.tsv \
    data/grch38_cpg_methylation.tsv \
    data/gnomad_nc_mutation_rates.tsv \
    data/vep_cds_all_possible_snvs.vcf \
    outputs/gnomad_pass_variants/all_pass_snvs.txt \
    outputs/nmd_annotations.tsv

# %% [markdown]
# ## Load datasets

# %%
# Define VCF headers and datatypes.
_header = ["chr", "pos", "id", "ref", "alt", "qual", "filter", "info"]

datatypes = defaultdict(lambda: "str")
datatypes.update({"pos": np.int32, "ac": np.int32, "an": np.int32})

# %%
# Retreive observed variants
obs = pd.read_csv(
    "../data/all_pass_snvs.txt",
    sep="\t",
    header=None,
    names=_header + ["ac", "an"],
    usecols=["chr", "pos", "ref", "alt", "ac", "an"],
    dtype=datatypes,
).assign(obs=1)

# %%
# Retreive VEP annotations of all possible SNVs
vep = pd.read_csv(
    "../data/vep_cds_all_possible_snvs.vcf",
    sep="\t",
    comment="#",
    header=None,
    names=_header,
    dtype=datatypes,
    usecols=["chr", "pos", "ref", "alt", "info"],
)

# %%
# Get enst
vep["enst"] = pd.Series([x.split("|", 3)[2] for x in vep["info"]])

# %%
# Get csq
syn = pd.Series(["synonymous" in x for x in vep["info"]])
mis = pd.Series(["missense" in x for x in vep["info"]])
non = pd.Series(["stop_gained" in x for x in vep["info"]])

vep.loc[syn, "csq"] = "synonymous"
vep.loc[mis, "csq"] = "missense"
vep.loc[non, "csq"] = "nonsense"

vep = vep.drop("info", axis=1).dropna()  # Keep only syn/mis/non variants

# %%
# Trinucleotide contexts
tri = pd.read_csv("../data/cds_trinucleotide_contexts.tsv", sep="\t", dtype=datatypes)

# %%
# gnomAD methylation data
meth = pd.read_csv(
    "../data/grch38_cpg_methylation.tsv",
    sep="\t",
    header=0,
    names=["ix", "chr", "pos", "alleles", "lvl"],
    usecols=["chr", "pos", "lvl"],
)

# %%
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
        "transversion":"non-CpG",
        "non-CpG transition":"non-CpG",
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

# %%
nmd = pd.read_csv(
    "../data/nmd_annotations.tsv",
    sep="\t",
    usecols=["chr", "pos", "transcript_id", "nmd_definitive"],
).rename(columns={"transcript_id": "enst", "nmd_definitive": "nmd"})

# %% [markdown]
# ## Merge annotations

# %%
# Merge VEP, context, NMD, and observed variant annotations
df = vep.merge(tri, how="left")
df = df.merge(nmd, how="left")
df = df.merge(obs, how="left").fillna(0)

# %%
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

# %% [markdown]
# ## Summarise the data

# %% [markdown]
# ### Rare synonymous variants for expectation model

# %%
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

# %% [markdown]
# ### Observed variants by region

# %%
# Transcripts
dfg = (
    df.groupby(["enst","csq", "variant_type"])
    .agg(
        n_pos=("pos", "count"),
        n_obs=("obs","sum"),
        mu=("mu","mean")
    )
    .reset_index()
)

# Save to output
dfg.to_csv("../outputs/observed_variants_stats_transcript.tsv", sep="\t", index=False)

# %%
# NMD regions
dfg = (
    df.groupby(["enst","csq", "variant_type", "nmd"])
    .agg(
        n_pos=("pos", "count"),
        n_obs=("obs","sum"),
        mu=("mu","mean")
    )
    .reset_index()
)

# Save to output
dfg.to_csv("../outputs/observed_variants_stats_nmd.tsv", sep="\t", index=False)

# %%
%%bash
# Upload outputs to UKB RAP
dx upload --destination outputs/ ../outputs/observed_variants_stats_transcript.tsv ../outputs/observed_variants_stats_nmd.tsv


