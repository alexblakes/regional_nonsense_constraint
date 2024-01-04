"""Module docstring."""

# Imports
from pathlib import Path

import pandas as pd

from src import constants as C
from src import setup_logger


# Logging
logger = setup_logger(Path(__file__).stem)


# Module constants
_NAMES = "location ref alt csq enst clinvar symbol".split()
_LOCATION_FIELDS = ["chr", "pos"]
_CLINVAR_FIELDS = "_index clinvar_symbol acmg review".split()
_CONSEQUENCES = [
    "missense_variant",
    "synonymous_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "inframe_deletion",
    "inframe_insertion",
    "stop_retained_variant",
    "start_retained_variant",
]  # The order matters
_COLUMNS = "chr pos ref alt symbol enst csq acmg review".split()


# Functions
def read_clinvar_vep(path):
    """Read VEP-annotated ClinVar variants."""

    logger.info(f"Reading VEP output.")

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=_NAMES,
        # nrows=1000000,  #! Testing
    )

    logger.info(f"VEP annotations: {len(df)}")

    return df


def split_fields(df):
    """Split fields in VEP output."""

    logger.info(f"Splitting VEP output fields.")

    df[_LOCATION_FIELDS] = df.location.str.split(":", expand=True)
    df[_CLINVAR_FIELDS] = df.clinvar.str.split("|", expand=True)

    df = df.drop(["location", "clinvar", "_index"], axis=1)

    logger.info(f"NaN values:\n{df.isna().sum()}")

    return df


def match_clinvar_symbol_to_vep_symbol(df):
    """Drop entries where ClinVar symbol does not match the VEP gene symbol."""

    df = df.query("symbol == clinvar_symbol").drop("clinvar_symbol", axis=1)

    logger.info(f"Entries after dropping unmatched gene symbols: {len(df)}")

    return df


def read_transcript_ids(path):
    """Get canonical transcript IDs."""

    return pd.read_csv(path, sep="\t", names="ensg enst symbol".split())


def filter_canonical_transcripts(df):
    """Filter for consequences in canonical transcripts."""

    canonical = read_transcript_ids(C.CANONICAL_CDS_GENE_IDS)["enst"].unique()

    df = df[df.enst.isin(canonical)]

    logger.info(f"Variants in canonical transcripts: {len(df)}")
    logger.info(
        f"Duplicates by variant and transcript ID: "
        f"{df.duplicated('chr pos ref alt enst'.split()).sum()}"
    )

    return df


def sanitise_vep_consequences(df):
    """Get one unambiguous consequence per variant."""

    # The order matters.
    for c in _CONSEQUENCES:
        df.loc[df.csq.str.contains(c), "csq"] = c

    df = df[df.csq.isin(_CONSEQUENCES)]

    _nl = "\n"  # Work-around as backslashes not allowed in f-strings
    logger.info(f"Consequences kept:\n{_nl.join(_CONSEQUENCES)}")
    logger.info(f"Variants after sanitising consequences: {len(df)}")
    logger.info(f"Consequence value counts:\n{df.csq.value_counts()}")

    return df


def main():
    """Run as script."""

    df = (
        read_clinvar_vep(C.CLINVAR_VEP)
        .pipe(split_fields)
        .pipe(match_clinvar_symbol_to_vep_symbol)
        .pipe(filter_canonical_transcripts)
        .pipe(sanitise_vep_consequences)
    )[_COLUMNS]

    return df  #! Testing


if __name__ == "__main__":
    df = main()
