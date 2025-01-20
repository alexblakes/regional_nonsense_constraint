"""Boilerplate code for most modules."""



import pandas as pd

import src

CANONICAL_TRANSCRIPTS = "data/interim/transcript_ids.tsv"
REGIONAL_NONSENSE_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"

logger = src.logger


def main():
    """Run as script."""

    canonical_transcripts = pd.read_csv(CANONICAL_TRANSCRIPTS).squeeze()

    rnc = (
        pd.read_csv(REGIONAL_NONSENSE_CONSTRAINT, sep="\t")
        .loc[lambda x: x["region"] != "transcript"]
        .loc[lambda x: x["constraint"] == "constrained"]
        .loc[lambda x: x["enst"].isin(canonical_transcripts)]
        .loc[:, ["enst", "pli", "loeuf"]]
        .drop_duplicates()
    )

    gnomad_weak = rnc.loc[lambda x: (x["pli"] <= 0.9) & (x["loeuf"] >= 0.6)]
    gnomad_absent = rnc.loc[lambda x: (x["pli"].isna()) & (x["loeuf"].isna())]

    logger.info(f"Regional constraint but weak gnomAD scores: {len(gnomad_weak)}")
    logger.info(f"Regional constraint but absent gnomAD scores: {len(gnomad_absent)}")

    return None


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
