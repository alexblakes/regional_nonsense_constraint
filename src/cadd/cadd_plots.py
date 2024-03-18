"""Plot CADD score summary data."""

# Imports
import logging
from pathlib import Path

import pandas as pd
import numpy as np

import src
from src import constants as C
from src import visualisation as vis
from src.visualisation import phylop_plots

# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_PALETTE = vis.color_palette("maps")
_METRIC = "CADD Phred"
_CSQS = ["Synonymous", "Missense", "Nonsense"]
_CONSTRAINT = ["constrained", "unconstrained"]

# Logging
logger = logging.getLogger(__name__)


# Functions
def read_data(path):
    return pd.read_csv(path, sep="\t")


def main():
    
    # Read data
    syn = read_data(C.STATS_CADD_SYN)
    mis = read_data(C.STATS_CADD_MIS)
    non = read_data(C.STATS_CADD_NON)
    
    return syn, mis, non


if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)
    main()
