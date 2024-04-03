"""Module docstring."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

import src
from src import constants as C
from src import visualisation as vis

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

logger = logging.getLogger(__name__)
plt.style.use(C.STYLE_DEFAULT)

df = pd.read_csv(C.STATS_GENE_SET_ENRICHMENT, sep="\t")
g = df.groupby(["background", "source", "query"])
g1 = g.get_group(("All genes", "HP", "gnomAD"))

fig, axs = plt.subplots(
    5, 3, figsize=(18 * C.CM, 30 * C.CM), sharex="col", layout="constrained"
)
axs = axs.flatten()

b = axs[0].barh(
    y=g1["enrichment_rank"],
    width=g1["enrichment"],
    label=g1["enrichment"],
    tick_label=g1["name"],
)
# axs[0].set_yticks(ticks=list(range(1, 11)), labels=g1["name"])
axs[0].bar_label(b, fmt="{:.2f}", padding=3)
axs[0].invert_yaxis()

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
