"""Module docstring."""

import itertools
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import src
from src import constants as C
from src import visualisation as vis

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

logger = logging.getLogger(__name__)

plt.style.use(C.STYLE_DEFAULT)
plt.style.use(C.COLOR_REGIONS)

df = pd.read_csv(C.STATS_GENE_SET_ENRICHMENT, sep="\t")
g = df.groupby(["background", "query", "source"])

fig, axs = plt.subplots(
    5, 3, figsize=(50 * C.CM, 30 * C.CM), layout="constrained"
)
axs = axs.flatten()

queries = ["gnomAD", "NMD target", "Start proximal", "Long exon", "Distal"]
sources = ["HP","GO:MF","GO:BP"]

for (q, s), ax, c in zip(itertools.product(queries, sources), axs, sns.color_palette()):
    
    data = g.get_group(("All genes", q, s)).copy()
        
    b = ax.barh(
        y=data["enrichment_rank"],
        width=data["enrichment"],
        tick_label=data["name"],
    )
    ax.bar_label(b, fmt="{:.2f}", padding=3)
    ax.invert_yaxis()

# Add HPO accession numbers to tick labels
# Add p-values to bar labels
# Colour the plot

plt.savefig("data/plots/gse_vs_all.png", dpi=600)
plt.close()


def main():
    """Run as script."""
    pass


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
