# %% [markdown]
# # Extract phyloP scores
# Get phyloP scores for all sites in canonical CDS.
# Use Zoonomia scores for 241 mammalian species.

# %% [markdown]
# ## Preliminaries
# Download data from UKB RAP storage

# %%
%%bash
dx download -o ../data/ data/241-mammalian-2020v2.bigWig
dx download -o ../outputs/ outputs/gencode_v39_canonical_cds_chr.bed

# %% [markdown]
# Install pyBigWig to conda environment

# %%
%%bash
conda install pybigwig -c conda-forge -c bioconda -y

# %%
# Import packages
import pyBigWig
import numpy as np
import pandas as pd

# %% [markdown]
# ## Script

# %%
# Read BigWig file of phyloP scores
bw = pyBigWig.open("../data/241-mammalian-2020v2.bigWig")
bw.header()

# %%
# Read bed file of canonical CDS regions
bed = pd.read_csv(
    "../outputs/gencode_v39_canonical_cds_chr.bed",
    sep="\t",
    header=None,
    names=["chr", "start", "end", "id", "score", "strand"],
    usecols=["chr", "start", "end"],
)

# Exlucde mitochondrial regions
bed = bed[bed["chr"] != "chrM"]

bed.head(3)

# %%
%%time
# Get phyloP annotations for each site in a CDS
## Extract annotations
phylop = bed.apply(lambda x: bw.intervals(x["chr"], x["start"], x["end"]), axis=1)

## Reformat the data
phylop.index = bed["chr"]
phylop = phylop.explode().dropna()
phylop = pd.DataFrame(
    [[*a] for a in phylop.values],
    columns=["start", "end", "phylop"],
    index=phylop.index,
)
phylop = (
    phylop.reset_index(drop=False).drop("start", axis=1).rename(columns={"end": "pos"})
).drop_duplicates()

# Print summary statistics
print(f"{len(phylop)} sites successfully annotated with phyloP scores.")

# %%
phylop.to_csv("../outputs/phylop_all_sites.tsv", sep="\t", index=False)

# %%
%%bash
dx upload --destination outputs/ ../outputs/phylop_all_sites.tsv

# %%



