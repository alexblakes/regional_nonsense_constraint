"""Parse the ClinVar summary text file.

Save tidied data as TSV and VCF for downstream annotation with VEP.
"""
from collections import defaultdict, OrderedDict

import pandas as pd

import src
from src import constants as C

logger = src.logger

FILE_IN = "data/raw/variant_summary.txt"
FILE_OUT = "data/interim/clinvar_parsed.tsv"
COLUMNS = OrderedDict(
    Name="name",
    ClinicalSignificance="acmg",
    OriginSimple="origin_simple",
    Assembly="assembly",
    Chromosome="chr",
    PositionVCF="pos",
    ReferenceAlleleVCF="ref",
    AlternateAlleleVCF="alt",
    ReviewStatus="review",
)
KEEP_CONTIGS = ["chr" + str(x) for x in list(range(1, 23))]
NULL_REVIEW = [
    "no assertion",
    "no interpretation",
]
NULL_ACMG = [
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
REPLACE_REVIEW = {
    "criteria provided, single submitter": "single",
    "criteria provided, multiple submitters, no conflicts": "multiple",
    "reviewed by expert panel": "expert",
    "practice guideline": "guideline",
}
REPLACE_ACMG = {
    "Uncertain significance": "VUS",
    "Likely benign": "B/LB",
    "Benign": "B/LB",
    "Benign/Likely benign": "B/LB",
    "Likely pathogenic": "P/LP",
    "Pathogenic": "P/LP",
    "Pathogenic/Likely pathogenic": "P/LP",
}


def read_data(path=FILE_IN):
    """Read ClinVar summary text file."""

    logger.info("Reading ClinVar variants.")

    dtypes = defaultdict(
        lambda: "category", Name="str", PositionVCF="int32", NumberSubmitters="int32"
    )

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=COLUMNS.keys(),
        low_memory=False,
        dtype=dtypes,
        na_values="na",
    ).rename(columns=COLUMNS)

    logger.info(f"ClinVar entries: {len(df)}\n\n" f"NaN values:\n{df.isna().sum()}")

    return df


@src.log_step
def filter_grch38(df):
    return df.query("assembly == 'GRCh38'").drop("assembly", axis=1)


def get_refseq_transcript_ids(df):
    df = df.assign(
        refseq_transcript=lambda x: x.name.str.extract(r"(^NM_[0-9]+)\.[0-9]+\(")
    ).drop("name", axis=1)

    logger.info(
        f"Entries without RefSeq transcript ID: {df.refseq_transcript.isna().sum()}"
    )

    return df


def extract_refseq_transcript_id(df):
    df["refseq_transcript"] = df.name.str.extract(r"(^NM_[0-9]+\.[0-9]+)\(")
    df = df.drop("name", axis=1)
    logger.info(f"Entries without{df.refseq_transcript.isna().sum()}")
    return df


def drop_nans(df):
    logger.info(f"NaN counts: {df.isna().sum()}")
    df = df.dropna()
    logger.info(f"Data shape: {df.shape}")
    return df


@src.log_step
def filter_contigs(df):
    return df.assign(
        chr=lambda x: x.chr.cat.rename_categories(lambda y: "chr" + str(y))
    ).loc[lambda z: z.chr.isin(KEEP_CONTIGS)]


@src.log_step
def filter_germline_origin(df):
    return df.loc[lambda x: x.origin_simple.str.contains("germline")].drop(
        "origin_simple", axis=1
    )


@src.log_step(display_args=True)
def filter_null_annotation(df, column: str, filter_terms: list):
    filter_string = "|".join(filter_terms)
    mask = df[column].str.lower().str.contains(filter_string)

    return df.loc[~mask]


def simplify_annotations(df, column: str, replace_dict: dict):
    df = df.replace({column: replace_dict})
    df[column] = df[column].cat.remove_unused_categories()

    logger.info(f"{column} value counts:\n{df[column].value_counts()}")

    return df


def drop_conflicting_acmg_annotations(df):
    """Drop variants with conflicting ACMG entries."""

    df = (
        df.drop_duplicates()
        .drop_duplicates("chr pos ref alt refseq_transcript acmg".split())
        .drop_duplicates("chr pos ref alt refseq_transcript".split(), keep=False)
    )

    logger.info(
        f"Data shape: {len(df)}\n\n"
        f"Duplicated by chr/pos/ref/alt: {df.duplicated('chr pos ref alt'.split()).sum()}\n\n"
        f"ACMG value counts:\n{df.acmg.value_counts()}\n\n"
        f"Review status value counts:\n{df.review.value_counts()}"
    )

    return df


def format_to_tsv(df):
    # Sorting required for VEP
    return (
        df.sort_values(["chr", "pos"])
        .reset_index()
        .loc[:, ["chr", "pos", "ref", "alt", "refseq_transcript", "acmg", "review"]]
    )


def main():
    """Run as script."""

    return (
        read_data(C.CLINVAR_VARIANT_SUMMARY)
        .pipe(filter_grch38)
        .pipe(get_refseq_transcript_ids)
        .pipe(drop_nans)
        .pipe(filter_contigs)
        .pipe(filter_germline_origin)
        .pipe(filter_null_annotation, "review", NULL_REVIEW)
        .pipe(simplify_annotations, "review", REPLACE_REVIEW)
        .pipe(filter_null_annotation, "acmg", NULL_ACMG)
        .pipe(simplify_annotations, "acmg", REPLACE_ACMG)
        .pipe(drop_conflicting_acmg_annotations)
        .pipe(format_to_tsv)
        .pipe(src.write_out, FILE_OUT)
    )


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()

# todo [X] get transcript id (NM_)
# todo [X] change "hgnc" to "symbol"
# todo [X] don't drop conflicting LP/P annotations, for example
# todo [X] check acmg simple
# todo [X] check origin simple
# todo [X] check n_submitters
