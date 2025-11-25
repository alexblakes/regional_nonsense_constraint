# %% [markdown]
# # Conservation
# Exploratory analysis of phyloP scores per region. 

# %%
import matplotlib.pyplot as plt
from matplotlib import ticker
import pandas as pd
from scipy import stats
import seaborn as sns

import src
from src import statistics_for_plots as sfp
from src import constants as C

# %% [markdown]
# ## Data munging

# %%
# Read data
df = pd.read_csv("data/interim/phylop_stats_per_region.tsv", sep="\t").pipe(
    sfp.sort_column, labels=C.NMD_REGIONS_DICT
)[
    [
        "symbol",
        "enst",
        "region",
        "constraint",
        "oe_ci_hi",
        "phylop_count",
        "phylop_median",
    ]
]

print(f"Data shape: {df.shape}")
df.head(3)

# %%
# Drop NaNs
df = df.dropna()
print(f"Data shape after dropping NaNs: {df.shape}")

# %%
# Drop genes symbols with multiple ENST IDs
print(f"Unique gene symbols: {df.symbol.nunique()}")
print(f"Unique ensts: {df.enst.nunique()}\n")

_uniq_symbol_enst = df[["symbol", "enst"]].drop_duplicates()
_symbol_dups = _uniq_symbol_enst.duplicated("symbol").sum()

print(
    f"There are {_symbol_dups} symbols with multiple transcripts.",
    "They are dropped.\n",
)

_kept_symbols = _uniq_symbol_enst["symbol"].drop_duplicates(keep=False)
df = df[df["symbol"].isin(_kept_symbols)].set_index("symbol")

print(f"Remaining unique gene symbols: {df.index.nunique()}")
print(f"Remaining unique ENSTs: {df.enst.nunique()}")

# %%
# Assign percentiles
df["oe_ci_hi_pct"] = df.groupby(["region"])["oe_ci_hi"].rank(method="min", ascending=False, pct=True)
df["phylop_median_pct"] = df.groupby(["region"])["phylop_median"].rank(method="min", ascending=True, pct=True)
df["pct_delta"] = df["oe_ci_hi_pct"] - df["phylop_median_pct"]

df.sample(5).sort_values(["region","oe_ci_hi"])

# %% [markdown]
# ## Find weakly conserved, strongly constrained regions

# %% [markdown]
# Filter for constrained regions with weak phyloP scores.

# %%
m1 = df["phylop_median_pct"] < 0.1
m2 = df["oe_ci_hi_pct"] > 0.9
m3 = df["constraint"] == "constrained"

lo_hi = df[m1 & m2 & m3]

print(
    f"Weakly conserved and strongly constrained transcript counts per region:\n"
    f"{lo_hi.groupby('region').size()}"
)

# %%
# Write out
lo_hi.pipe(src.write_out, "data/statistics/constraint_conservation_discrepant_genes.tsv", index=True)

# %%
lo_hi.sort_values(["region","oe_ci_hi_pct"], ascending=False)

# %% [markdown]
# Many of these regions are poorly covered, and therefore have strong oe_ci_95 scores, but have an "indeterminate" constraint annotation.
# It will be sensible to filter for constrained regions only.
# 
# One example is the last exon of NANOG (ENST00000229307, OE95=0.43, mean phyloP=0.2).
# This appears to be a mammalian expansion of the CDS. 
# The highly conserved homeodomains occur in two smaller internal exons.
# A quick scan of the literature suggests that the C-terminal domain of human NANOG contains two highly potent transactivating domains.
# Although these lack any homology or structural resemblance to known domains and are poorly modelled by alphafold.
# 
# ? Immune genes e.g. CFH (ENST00000367429), PTPRC (ENST00000442510), 
# 
# Non-conserved long-exon region of ICE1
# 
# The central exons of TCOF (Treacher-Collins Syndrome) are poorly conserved.
# 
# The start-proximal region of USP17L2 is strongly constrained but weakly conserved.
# 
# Four ZNF proteins have strong distal constraint but weak conservation.
# (But they don't look particularly interesting.)
# 
# **Jamie flags CFH as a flagship gene at the interface between rare and common eye disorders.**

# %% [markdown]
# ## Find the strongest outliers in each region

# %%
constrained = df[df["constraint"] == "constrained"]
constrained.groupby("region")["pct_delta"].nlargest(5)

# %% [markdown]
# ARID3B looks an interesting case. It has a long first exon which is lowly conserved, but totally depleted for nonsense variants in gnomAD.
# The remainder of the CDS is more highly conserved.
# 
# Similarly RNF6. The first two coding exons have regions of weak conservation, but are totally depleted for nonsense or frameshift variants.
# 
# 

# %% [markdown]
# ## Distribution of centiles
# It will be interesting to look at the deviation between constraint and conservation within each region.
# For example, distal regions seem enriched among regions with the strongest constraint but weakest conservation (above).
# Is this a general property, or just something we see in a handful of examples?

# %%
plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])


def standard_figure(**kwargs):
    kwargs.setdefault("nrows", 1)
    kwargs.setdefault("ncols", 4)
    kwargs.setdefault("figsize", (18 * C.CM, 4.5 * C.CM))
    kwargs.setdefault("layout", "constrained")

    return plt.subplots(**kwargs)


fig, axs = standard_figure()

for ax, (name, group), color in zip(axs, df.groupby("region"), sns.color_palette()[1:]):
    # Look at phyloP medians here...
    statistic = group["oe_ci_hi_pct"] - group["phylop_median_pct"]

    sns.kdeplot(ax=ax, x=statistic, color=color, fill=True, cut=0)

    ax.set_xlabel("Delta percentile")
    ax.set_title(f"{name}\nn={len(group)}", pad=20)
    ax.axvline(x=statistic.mean(), color="black", ls="--")

    # Annotations for clarity
    ax.annotate(
        "",
        xy=(0.05, 0.9),
        xytext=(0.45, 0.9),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="->"),
    )
    ax.annotate(
        "More\nconserved",
        xy=(0.25, 1),
        xytext=(0.25, 1),
        textcoords="axes fraction",
        ha="center",
        ma="center",
    )

    ax.annotate(
        "",
        xy=(0.95, 0.9),
        xytext=(0.55, 0.9),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="->"),
    )
    ax.annotate(
        "More\nconstrained",
        xy=(0.75, 1),
        xytext=(0.75, 1),
        textcoords="axes fraction",
        ha="center",
        ma="center",
    )

# %% [markdown]
# There is no glaring difference between the regions

# %%
fig, axs = standard_figure(sharex=True)

for ax, (name, group), color in zip(axs, df.groupby("region"), sns.color_palette()[1:]):
    # Limit to constrained regions
    group = group[group["constraint"] == "constrained"]

    # Look at phyloP medians here...
    statistic = group["oe_ci_hi_pct"] - group["phylop_median_pct"]

    print(
        f"{name} summary statistics:\n{stats.describe(statistic, nan_policy='omit')}\n"
    )
    sns.kdeplot(ax=ax, x=statistic, color=color, fill=True, cut=0)

    ax.annotate(
        "",
        xy=(1, 0.7),
        xytext=(0.55, 0.7),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="->"),
    )
    ax.annotate("More\nconstrained", xy=(1, 0.75), xytext=(1,0.75), textcoords="axes fraction", ha="right", ma="left")

    ax.set_xlabel("Delta percentile")
    ax.set_title(f"{name} (constrained)\nn={len(group):,}")

    ax.axvline(statistic.mean(), color="black", ls="--")
    ax.annotate(f"{statistic.mean():.2}", xy=(statistic.mean(),0), xytext=(3, 3), textcoords="offset points")

plt.savefig("data/plots/conservation/constraint_conservation_pct_deltas.png", dpi=600)
plt.savefig("data/plots/conservation/constraint_conservation_pct_deltas.svg")

# %% [markdown]
# Constrained long exon and distal regions have higher mean delta percentile.
# These plots could be expressed as box and whiskers.
# Significance could be test with a T test.
# 
# The high level trend appears to be that constrained regions have positive percentile differences.
# But this is an artefact of selecting for constrained regions only. 
# In these regions, there is only one direction in which conserved values can go if they are to be different from constrained values.

# %% [markdown]
# ## Comparing raw scores

# %% [markdown]
# First, create scatter plots of constraint vs phyloP raw scores.

# %%
def hexbin_plot(x, y, ax=None, **kwargs):
    kwargs.setdefault("xscale", "log")
    kwargs.setdefault("bins", "log")
    kwargs.setdefault("cmap", "Blues_r")
    # kwargs.setdefault("edgecolors", "none")
    kwargs.setdefault("linewidths", 0)
    
    ax = ax or plt.gca()

    hb = ax.hexbin(
        x,
        y,
        **kwargs,
    )

    return hb


def customise_hexbin(x, y, title, ax=None):
    ax = ax or plt.gca()

    ax.set_xlabel("OE95")
    ax.set_ylabel("Median phyloP")
    ax.set_title(f"{title}\nn={len(x):,}", pad=15)

    rho = stats.spearmanr(x, y).statistic
    ax.annotate(
        f"rho={rho:.2f}",
        xy=(1, 1),
        xycoords=("axes fraction"),
        ha="right",
        va="bottom",
        ma="right",
    )

    return ax


fig, axs = standard_figure(figsize=(18 * C.CM, 5 * C.CM), sharex=True, sharey=True)

for ax, (name, group) in zip(axs, df.groupby("region")):
    x = group["oe_ci_hi"]
    y = group["phylop_median"]

    hb = hexbin_plot(x, y, ax)
    customise_hexbin(x, y, name, ax)

    # Simplify the tick marks
    ax.xaxis.set_major_formatter('{:g}'.format)

cb = fig.colorbar(hb, ax=axs[3], label="Regions", format=ticker.ScalarFormatter())

# %% [markdown]
# Constraint in long exon and distal regions is less strongly correlated with conservation.
# 
# What about in constrained regions only?

# %%
fig, axs = standard_figure(figsize=(18 * C.CM, 5 * C.CM), sharex=True, sharey=True)

for ax, (name, group) in zip(axs, df.groupby("region")):
    # Subset to constrained regions only
    group = group[group.constraint == "constrained"]
    
    x = group["oe_ci_hi"]
    y = group["phylop_median"]

    hb = hexbin_plot(x, y, ax)
    customise_hexbin(x, y, name, ax)

    # Simplify the tick marks
    ax.xaxis.set_major_locator(ticker.FixedLocator([0.1, 0.3, 0.6]))
    ax.xaxis.set_major_formatter("{:g}".format)

cb = fig.colorbar(hb, ax=axs[3], label="Regions", format=ticker.ScalarFormatter())

# %% [markdown]
# The story is less clearcut, but again distal regions correlate most poorly.
# 
# Perhaps comparing percentiles in constrained regions will be clearer?

# %% [markdown]
# ## Comparing percentile scores 

# %%
fig, axs = standard_figure(figsize=(18 * C.CM, 5 * C.CM), sharex=True, sharey=True)

for ax, (name, group) in zip(axs, df.groupby("region")):
   
    x = group["oe_ci_hi_pct"]
    y = group["phylop_median_pct"]

    hb = hexbin_plot(x, y, ax, xscale="linear", bins="log")
    customise_hexbin(x, y, name, ax)

cb = fig.colorbar(hb, ax=axs[3], label="Regions", format=ticker.ScalarFormatter())

# %% [markdown]
# No clear trend emerges, because the rank normalised scores are uniformly distributed.

# %%
fig, axs = standard_figure(figsize=(18 * C.CM, 5 * C.CM), sharey=True)

for ax, (name, group) in zip(axs, df.groupby("region")):
    # Subset to constrained regions only
    group = group[group.constraint == "constrained"]
       
    x = group["oe_ci_hi_pct"]
    y = group["phylop_median_pct"]

    hb = hexbin_plot(x, y, ax, xscale="linear", bins="log")
    customise_hexbin(x, y, name, ax)

cb = fig.colorbar(hb, ax=axs[3], label="Regions", format=ticker.ScalarFormatter())

# %% [markdown]
# It's also not a particularly clear message after subsetting to constrained genes.

# %% [markdown]
# ## Normalisation to find outliers

# %% [markdown]
# ## Conclusion
# Of the genes prioritised with this approach, NANOG looks the most interesting.
# There is no regional pattern for divergence between constraint and phylop scores on aggregate.
# But for constrained regions, long exon and distal regions tend to have greater disparity between constraint and conservation. 

# %% [markdown]
# ## To do
# - [X] Distributions of constraint  / median conservation percentile deltas, per region
# - [ ] Highlight outliers with the smallest "sum of percentiles" - maybe 5-10 genes in each region.

# %% [markdown]
# 


