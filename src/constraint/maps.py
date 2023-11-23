
# # Observed variants
# This script identifies the observed and possible variants in UKB. 
# 
# Variants are grouped by variant context, transcript, or NMD-region. The script aggregates the number observed, the number possible, and the mean mutability for each grouping. The summary data are saved to a .tsv outputs.


# ## Preliminaries


# ### Import modules


import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as _stats
import statsmodels.formula.api as smf
from statsmodels.stats.proportion import proportions_ztest

sns.set_context("talk")


# ### Download datasets from UKB RAP

dx download \
    -f \
    -o ../data/ \
    data/cds_trinucleotide_contexts.tsv \
    data/grch38_cpg_methylation.tsv \
    data/gnomad_nc_mutation_rates.tsv \
    data/vep_cds_all_possible_snvs.vcf \
    outputs/gnomad_pass_variants/all_pass_snvs.txt \
    outputs/nmd_annotations.tsv


# ## Load datasets


# Define VCF headers and datatypes.
_header = ["chr", "pos", "id", "ref", "alt", "qual", "filter", "info"]

datatypes = defaultdict(lambda: "str")
datatypes.update({"pos": np.int32, "ac": np.int32, "an": np.int32})


# Retreive observed variants
obs = pd.read_csv(
    "../data/all_pass_snvs.txt",
    sep="\t",
    header=None,
    names=_header + ["ac", "an"],
    usecols=["chr", "pos", "ref", "alt", "ac", "an"],
    dtype=datatypes,
).assign(obs=1)


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


# Get enst
vep["enst"] = pd.Series([x.split("|", 3)[2] for x in vep["info"]])


# Get csq
syn = pd.Series(["synonymous" in x for x in vep["info"]])
mis = pd.Series(["missense" in x for x in vep["info"]])
non = pd.Series(["stop_gained" in x for x in vep["info"]])

vep.loc[syn, "csq"] = "synonymous"
vep.loc[mis, "csq"] = "missense"
vep.loc[non, "csq"] = "nonsense"

vep = vep.drop("info", axis=1).dropna()  # Keep only syn/mis/non variants


# Trinucleotide contexts
tri = pd.read_csv("../data/cds_trinucleotide_contexts.tsv", sep="\t", dtype=datatypes)


# gnomAD methylation data
meth = pd.read_csv(
    "../data/grch38_cpg_methylation.tsv",
    sep="\t",
    header=0,
    names=["ix", "chr", "pos", "alleles", "lvl"],
    usecols=["chr", "pos", "lvl"],
)


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

# Mutation rates are only available for 32 codons. We need to reverse-complement for the remainder.
complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
# Replace ref and alt alleles
_mu = mu.copy().replace(complement)
# Reverse-complement trinucleotide contexts
_mu["tri"] = pd.Series(["".join([complement[y] for y in x])[::-1] for x in mu.tri])
mu = pd.concat([mu, _mu])


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


# ### Synonymous proportion singletons by variant context


def get_ps(dfg):
    """
    Get mean mutability and proportion of singletons.
    """
    mu = dfg["mu"].mean()
    ns = dfg["ac"].apply(lambda x: (x==1).sum()).rename("n_singletons")
    no = dfg["ac"].count().rename("n_obs")
    
    ps = pd.concat([mu, ns, no], axis=1).reset_index()
    
    ps["ps"] = ps["n_singletons"] / ps["n_obs"]
    
    return ps


# Subset to synonymous variants only
syn = df[df["csq"] == "synonymous"].copy()

# Mask contexts in which a synonymous variant is generally not possible.
# (NB synonymous variants in these contexts can only occur at exon-intron junctions)
m1 = (syn.tri == "AGT") & ((syn.alt == "C") | (syn.alt == "T"))
m2 = (syn.tri == "AAT") & ((syn.alt == "C") | (syn.alt == "T"))
m3 = (syn.tri == "ACT") & ((syn.alt == "G") | (syn.alt == "A"))
m4 = (syn.tri == "ATT") & ((syn.alt == "G") | (syn.alt == "A"))

# Mask variants which are not observed
m5 = syn["ac"] == 0

# Apply filters
syn_obs = syn[~(m1 | m2 | m3 | m4 | m5)]

# Get proportion of singletons per context for observed variants
syn_g = syn_obs.groupby(["tri", "ref", "alt", "variant_type", "lvl"])

ps = get_ps(syn_g)

ps.to_csv("../outputs/proportion_singletons_synonymous_by_context.tsv", sep="\t", index=False)


# ### Proportion of singletons by csq


m1 = df["ac"] != 0
m2 = df["variant_type"] != "CpG"
m3 = df["csq"] == "nonsense"


# #### Synonymous, missense, and nonsense


# CpG only
ps_csq_cpg = get_ps(df[m1 & ~m2].groupby("csq"))

# Exclude CpG
ps_csq_no_cpg = get_ps(df[m1 & m2].groupby("csq"))


# #### Nonsense, by NMD region


# CpG only
ps_region_cpg = get_ps(df[m1 & m3 & ~m2].groupby("nmd"))

# Exclude CpG
ps_region_no_cpg = get_ps(df[m1 & m2 & m3].groupby("nmd"))


# ### Combine CSQ and Region results


# Reformat column names
ps_region_cpg = ps_region_cpg.rename(columns={"nmd":"csq"})
ps_region_no_cpg = ps_region_no_cpg.rename(columns={"nmd":"csq"})

for a, b, c in zip([ps_csq_cpg, ps_csq_no_cpg], [ps_region_cpg, ps_region_no_cpg], ["cpg","no_cpg"]):
    ps = pd.concat([a,b])
    ps.to_csv(f"../outputs/proportion_singletons_by_csq_and_region_{c}.tsv", sep="\t", index=False)
    print(f"{c}\n{ps}\n\n")




