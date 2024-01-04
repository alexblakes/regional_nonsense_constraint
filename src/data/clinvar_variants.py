"""Parse the ClinVar summary text file.

Save tidied data as TSV and VCF for downstream annotation with VEP.
"""

# Imports
from pathlib import Path

import pandas as pd

from src import constants as C
from src import setup_logger
from src.data import clinvar_vep_tidy


# Logging
logger = setup_logger(Path(__file__).stem)


# Module constants
_USECOLS = [
    "Type",
    "GeneSymbol",
    "ClinicalSignificance",
    "Assembly",
    "Chromosome",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "ReviewStatus",
]

_NAMES = [
    "type",
    "hgnc",
    "acmg",
    "assembly",
    "chr",
    "pos",
    "ref",
    "alt",
    "review",
]

_CHROM = ["chr" + str(x) for x in list(range(1, 23))] + ["X", "Y"]

_NULL_REVIEW = [
    "no assertion",
    "no interpretation",
]

_NULL_ACMG = [
    "not provided",
    "drug response",
    "other",
    "risk",
    "low penetrance",
    "conflicting",
    "affects",
    "association",
    "protective",
    "confers sensitivity",
]

_ACMG = [
    "uncertain significance",
    "likely benign",
    "benign",
    "likely pathogenic",
    "pathogenic",
]


# Functions
def read_clinvar_summary(path):
    """Read ClinVar summary text file."""

    logger.info("Reading ClinVar variants.")

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=_USECOLS,
        low_memory=False,
        dtype={"Chromosome": str},
        na_values="na",
    ).rename(columns={x: y for x, y in zip(_USECOLS, _NAMES)})

    logger.info(f"ClinVar entries: {len(df)}")
    logger.info(f"NaN values:\n{df.isna().sum()}")

    # Drop entries with missing data.
    df = df.dropna()

    logger.info(f"Entries after dropping NaNs: {len(df)}")
    logger.info(
        f"Duplicated by variant: {df.duplicated('chr pos ref alt'.split()).sum()}"
    )

    return df


def filter_grch38(df):
    """Filter for entries in GRCh38."""

    df = df.query("assembly == 'GRCh38'").drop("assembly", axis=1)

    # Add "chr" prefix to chromosome numbers
    df["chr"] = "chr" + df["chr"]

    return df


def filter_major_contigs(df):
    """Filter for variants on the major chromosomes."""

    df = df.query(f"chr.isin({_CHROM})")

    logger.info(f"Variants on major contigs: {len(df)}")

    return df


def filter_null_annotations(df):
    """Exclude variants with irrelevant review or ACMG annotations."""

    # Create masks for filtering
    m1 = ~df.review.str.lower().str.contains("|".join(_NULL_REVIEW))
    m2 = ~df.acmg.str.lower().str.contains("|".join(_NULL_ACMG))

    # Filter for irrelevant review or ACMG values
    df = df[m1 & m2]

    logger.info(
        f"ClinVar variants with relevant review and ACMG annotations: {len(df)}"
    )
    logger.info(
        f"Duplicated by variant: {df.duplicated('chr pos ref alt'.split()).sum()}"
    )

    return df


def drop_conflicting_acmg_annotations(df):
    """Drop variants with conflicting ACMG entries."""

    df = (
        df.replace(
            {
                "Benign/Likely benign": "Likely benign",
                "Pathogenic/Likely pathogenic": "Likely pathogenic",
            }
        )
        .drop_duplicates()
        .drop_duplicates("chr pos ref alt hgnc acmg".split())
        .drop_duplicates("chr pos ref alt hgnc".split(), keep=False)
    )

    logger.info(
        f"Variants remaining after dropping conflicting ACMG annotations: {len(df)}"
    )
    logger.info(f"ACMG value counts:\n{df.acmg.value_counts()}")
    logger.info(f"Review status value counts:\n{df.review.value_counts()}")

    return df


def rationalise_gene_symbols(df):
    """Keep only one HGNC gene symbol per variant."""

    # Get HGNC symbols for canonical transcripts.
    symbols = clinvar_vep_tidy.read_transcript_ids(C.CANONICAL_CDS_GENE_IDS).symbol

    # Explode and filter any nested gene symbols.
    df["hgnc"] = df["hgnc"].str.split(";")
    df = df.explode("hgnc")
    df = df[df.hgnc.isin(symbols)]

    logger.info(f"ClinVar entries with HGNC symbols: {len(df)}")

    return df


def format_to_tsv(df):
    df = df.sort_values(["chr", "pos"]).reset_index()  # Required for VEP
    return df[["chr", "pos", "ref", "alt", "hgnc", "acmg", "review"]]


def format_to_vcf(df):
    """Bespoke function to reformat ClinVar data to VCF."""

    # Placeholder values for qual, filter, and info columns.
    df = df.assign(
        qual=".",
        _filter=".",
        _info=".",
    )

    # Put ClinVar annotations in the id column.
    df["_id"] = [
        "|".join([str(a), b, c, d])
        for a, b, c, d in zip(df.index, df.hgnc, df.acmg, df.review)
    ]

    # Select columns for VCF format.
    df = df["chr pos _id ref alt qual _filter _info".split()]

    return df


def write_tsv(df, path, **kwargs):
    logger.info(f"Writing to output at {path}")
    df.to_csv(path, sep="\t", index=False, **kwargs)
    pass


def main():
    """Run as script."""

    clinvar = (
        read_clinvar_summary(C.CLINVAR_VARIANT_SUMMARY)
        .pipe(filter_grch38)
        .pipe(filter_major_contigs)
        .pipe(filter_null_annotations)
        .pipe(drop_conflicting_acmg_annotations)
        .pipe(rationalise_gene_symbols)
        .pipe(format_to_tsv)
    )
    vcf = format_to_vcf(clinvar)

    write_tsv(clinvar, C.CLINVAR_SELECTED_TSV)
    write_tsv(vcf, C.CLINVAR_SELECTED_VCF, header=False)

    return clinvar


if __name__ == "__main__":
    clinvar = main()
