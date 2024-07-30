"""Plot CADD scores in constrained and unconstrained regions."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

import src
from src import constants as C
from src.visualisation import orthogonal_metrics as om

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/statistics/orthogonal_metrics_cadd.tsv.gz"
_PNG = "data/plots/orthogonal_metrics/cadd.png"
_SVG = "data/plots/orthogonal_metrics/cadd.svg"

logger = logging.getLogger(__name__)


def main():
    """Run as script."""

    plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])
    fig, ax = plt.subplots(1, 1, figsize=(6 * C.CM, 6 * C.CM), layout="constrained")

    cadd = om.read_data(_FILE_IN).pipe(
        om.sort_data, categories=[x for x in C.REGION_LABELS if x != "Full CDS"]
    )

    om.grouped_boxplot(cadd, ax, y="cadd_phred")
    om.customise_plot(ax, ylabel="CADD phred", legend=True)
    om.recolor(ax, sns.color_palette()[1:])

    logger.info("Saving plot.")
    plt.savefig(_PNG, dpi=600)
    plt.savefig(_SVG)
    plt.close("all")

    return ax


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
