# %%
# Plot the proportion of NMD regions which overlap with a protein domain.

import pandas as pd
import pandas_checks as pdc
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
from statannotations.Annotator import Annotator
from statsmodels.stats import proportion

import vis

FILE_IN = "data/statistics/proportion_nmd_regions_overlapping_pfam_domains.tsv"
PNG = "data/plots/pfam_domains/nmd_regions_overlapping_pfam_domains.png"
SVG = "data/plots/pfam_domains/nmd_regions_overlapping_pfam_domains.svg"

# %%
df = pd.read_csv(FILE_IN, sep="\t").check.function()

# %%
# Overlap of NMD regions with protein domains (ignoring constraint)
df.groupby("region").agg(count=("count","sum"), total=("total","sum")).assign(pro=lambda x: x["count"] / x["total"])

# %%
# Quantify statistical difference between constrained and unconstrained regions
dfg_region = df.groupby("region", sort=False)

print("Two-sample Z test of proportions, constrained vs unconstrained")

for name, data in dfg_region:
    count = data["count"]
    nobs = data["total"]
    z, p = proportion.proportions_ztest(count, nobs, alternative="two-sided")
    print(f"{name}:\n\tZ={z}\n\tP={p}")

# %%
fig, ax = plt.subplots(figsize=(8 / 2.54, 8 / 2.54), layout="constrained")

# Plot grouped bar chart
dark_colors = ["#B2182B", "#98CAE1", "#4A7BB7", "#293C83"]
light_colors = [vis.adjust_alpha(c, 0.35) for c in dark_colors]

vis.plot.grouped_vertical_bar(
    df,
    "region",
    "constraint",
    "prop",
    ax=ax,
    bar_label_fmt="{:.2f}",
    err_columns=["err_lo", "err_hi"],
    edgecolor="k",
    colors=[light_colors, dark_colors],
)
## Legend
legend = ax.legend()
for t in legend.get_texts():
    t.set_color("k")

## Axis labels
ax.set_ylabel("Proportion of regions\ncontaining a Pfam domain")
vis.utils.rotate_tick_labels("x", ax=ax)

# Add stats annotation
pairs = [
    (("NMD target", "Unconstrained"), ("NMD target", "Constrained")),
    (("Start proximal", "Unconstrained"), ("Start proximal", "Constrained")),
    (("Long exon", "Unconstrained"), ("Long exon", "Constrained")),
    (("Distal", "Unconstrained"), ("Distal", "Constrained")),
]
pvals = [
    1.11321727493653e-19,
    0.6539163322708069,
    0.1446197470023705,
    0.00016691760093622387,
]
annot = Annotator(
    ax=ax, pairs=pairs, data=df, x="region", y="prop", hue="constraint", plot="boxplot"
)
vis.utils.configure_annotator(
    annot, test=None, pvalue_format_string="{:.3G}"
).set_pvalues(pvals).annotate()

# Save figure
plt.savefig(PNG, dpi=600)
plt.savefig(SVG)


