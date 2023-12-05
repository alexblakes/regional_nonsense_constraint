"""Tidy pext scores prior to liftover."""

# Imports
from pathlib import Path

import pandas as pd

from src import constants as C
from src import setup_logger


# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions

def main():
    """Run as script."""

    df = pd.read_csv("C.", sep="\t").dropna().reset_index()
    

    return df #! Testing

df = pd.read_csv("../outputs/pext_37.tsv", sep="\t").dropna().reset_index()


_ = pd.DataFrame([x.split(":") for x in df.locus], columns=["chr","end"])

df["chr"] = "chr" + _["chr"]
df["end"] = _["end"].astype(int)
df["start"] = df["end"] - 1
df = df.rename(columns={"ensg":"name", "mean_proportion":"score"})

df = df[["chr","start","end","name","score"]]

df.to_csv("../outputs/pext_37.bed", sep="\t", index=False, header=False)


