"""Explore the relationship between nonsense O/E values and LOEUF scores."""



import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import src
from src import constants as C


_FILE_IN = "data/final/regional_nonsense_constraint.tsv"
OE_SVG = "data/plots/constraint/oe_vs_loeuf.svg"
OE_PNG = "data/plots/constraint/oe_vs_loeuf.png"
OE_CI_SVG = "data/plots/constraint/oe_ci_hi_vs_loeuf.svg"
OE_CI_PNG = "data/plots/constraint/oe_ci_hi_vs_loeuf.png"


logger = src.logger


def main():
    """Run as script."""

    df = pd.read_csv(_FILE_IN, sep="\t")

    return df


if __name__ == "__main__":
    src.add_log_handlers()
    main()
