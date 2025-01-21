"""Parse the ClinVar summary text file.

Save tidied data as TSV and VCF for downstream annotation with VEP.
"""
from collections import defaultdict

import pandas as pd

import src
from src import constants as C

logger = src.logger

CLINVAR_VARIANT_SUMMARY = "data/raw/variant_summary.txt"

_COLUMNS = {
    "Name": "name",
    "GeneSymbol": "hgnc",
    "ClinicalSignificance": "acmg",
    "ClinSigSimple": "acmg_simple",
    "OriginSimple": "origin_simple",
    "Assembly": "assembly",
    "Chromosome": "chr",
    "PositionVCF": "pos",
    "ReferenceAlleleVCF": "ref",
    "AlternateAlleleVCF": "alt",
    "ReviewStatus": "review",
    "NumberSubmitters": "n_submitters",
}
_CHROMS = [
    "chr" + str(x) for x in list(range(1, 23))
]  # Autosomes and major contigs only
_NULL_REVIEW = [
    "no assertion",
    "no interpretation",
]
_REVIEW = [
    "criteria provided, single submitter",
    "criteria provided, multiple submitters, no conflicts",
    "reviewed by expert panel",
    "practice guideline",
]
_REVIEW_BRIEF = "single multiple expert guideline".split()
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
    "Uncertain significance",
    "Likely benign",
    "Benign",
    "Likely pathogenic",
    "Pathogenic",
]
_ACMG_BRIEF = "VUS LB B LP P".split()


def read_data(path=CLINVAR_VARIANT_SUMMARY):
    """Read ClinVar summary text file."""

    logger.info("Reading ClinVar variants.")

    dtypes = defaultdict(lambda: "category")
    dtypes.update(name=str,pos=int, n_submitters=int)

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=_COLUMNS.keys(),
        low_memory=False,
        dtype=dtypes,
        na_values="na",
        nrows=100,
    ).rename(columns=_COLUMNS)

    logger.info(f"ClinVar entries: {len(df)}\n" f"NaN values:\n{df.isna().sum()}")

    return df

    # Drop entries with missing data.
    df = df.dropna()

    logger.info(f"Entries after dropping NaNs: {len(df)}")
    logger.info(
        f"Duplicated by chr/pos/ref/alt: {df.duplicated('chr pos ref alt'.split()).sum()}"
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

    df = df.query(f"chr.isin({_CHROMS})")

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
        f"Duplicated by chr/pos/ref/alt: {df.duplicated('chr pos ref alt'.split()).sum()}"
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
    logger.info(
        f"Duplicated by chr/pos/ref/alt: {df.duplicated('chr pos ref alt'.split()).sum()}"
    )
    logger.info(f"ACMG value counts:\n{df.acmg.value_counts()}")
    logger.info(f"Review status value counts:\n{df.review.value_counts()}")

    return df


def replace_annotations(df, col, replace_list, value_list):
    replace_dict = {a: b for a, b in zip(replace_list, value_list)}
    df = df.replace({col: replace_dict})
    logger.info(f'Value counts in "{col}":\n{df[col].value_counts()}')

    return df


def rationalise_gene_symbols(df):
    """Keep only one HGNC gene symbol per variant."""

    # Get HGNC symbols for canonical transcripts.
    symbols = pd.read_csv(C.CANONICAL_CDS_GENE_SYMBOLS, header=None).squeeze()

    # Explode and filter any nested gene symbols.
    df["hgnc"] = df["hgnc"].str.split(";")
    df = df.explode("hgnc")
    df = df[df.hgnc.isin(symbols)]

    logger.info(f"ClinVar entries with HGNC symbols: {len(df)}")

    return df


def format_to_tsv(df):
    df = df.sort_values(["chr", "pos"]).reset_index()  # Required for VEP
    return df[["chr", "pos", "ref", "alt", "hgnc", "acmg", "review"]]


def main():
    """Run as script."""

    clinvar = (
        read_data(C.CLINVAR_VARIANT_SUMMARY)
        # .pipe(filter_grch38)
        # .pipe(filter_major_contigs)
        # .pipe(filter_null_annotations)
        # .pipe(drop_conflicting_acmg_annotations)
        # .pipe(replace_annotations, "acmg", _ACMG, _ACMG_BRIEF)
        # .pipe(replace_annotations, "review", _REVIEW, _REVIEW_BRIEF)
        # .pipe(rationalise_gene_symbols)
        # .pipe(format_to_tsv)
    )

    # logger.info("Writing to output.")
    # clinvar.to_csv(C.CLINVAR_SELECTED_TSV, sep="\t", index=False, header=False)

    return clinvar


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
