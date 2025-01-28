"""Boilerplate code for most modules."""

import re

import pandas as pd

import src

FILE_IN = "data/raw/gwas_catalog_all_association_v1.0.2.tsv"
COLUMNS = {
    "DISEASE/TRAIT": "trait",
    "CHR_ID": "chrom",
    "CHR_POS": "pos",
    "STRONGEST SNP-RISK ALLELE": "strongest_allele",
    "SNPS": "snps",
    "CONTEXT": "csq",
    "RISK ALLELE FREQUENCY": "af",
    "P-VALUE": "p",
}
logger = src.logger


@src.log_step
def read_data(path=FILE_IN):
    df = pd.read_csv(path, sep="\t", usecols=COLUMNS, low_memory=False).rename(
        columns=COLUMNS
    )

    logger.info(f"Data shape: {df.shape}\n\n" f"Data info:\n{df.isna().sum()}")

    return df


@src.log_step
def drop_nans_in_chrom_pos(df):
    return df.dropna(subset=["chrom", "pos"])


@src.log_step
def drop_haplotypes(df):
    return df[~df.chrom.str.contains(";")]


@src.log_step
def drop_snp_interactions(df):
    df = df[~df.chrom.str.contains(r" x ")]
    logger.info(f"Data shape: {df.shape}\n\n" f"Unique chroms: {df.chrom.nunique()}")
    return df


@src.log_step
def drop_unknown_alleles(df):
    return df[~df.strongest_allele.str.contains(r"\?")]


@src.log_step
def drop_missing_strongest_alleles(df):
    return df[df.strongest_allele.str.contains(r"-")]


@src.log_step
def get_alt_allele(df):
    pattern = re.compile(r".*-(.*)")
    return df.assign(alt=lambda x: x.strongest_allele.str.extract(pattern)).drop(
        columns="strongest_allele"
    )


@src.log_step
def drop_non_ATCG_alt_alleles(df):
    pattern = re.compile(r"[^ATCG]")
    mask_non_atcg = df.alt.str.contains(pattern)

    logger.info(f"Non-ATCG value counts:\n{df[mask_non_atcg].alt.value_counts()}")

    return df[~mask_non_atcg]


def pos_to_int(df):
    return df.astype({"pos": "int32"})


def drop_interaction_csqs(df):
    # Two SNPs where csq == "intergenic_variant x intron_variant"
    df = df[~df.csq.fillna("").str.contains(r" x ")]
    logger.info(
        f"Data shape: {df.shape}\n\n"
        f"Consequence value counts:\n{df.csq.value_counts(dropna=False)}"
    )
    return df


def main():
    """Run as script."""
    return (
        read_data()
        .pipe(drop_nans_in_chrom_pos)
        .pipe(drop_haplotypes)
        .pipe(drop_snp_interactions)
        .pipe(drop_unknown_alleles)
        .pipe(drop_missing_strongest_alleles)
        .pipe(get_alt_allele)
        .pipe(drop_non_ATCG_alt_alleles)
        .pipe(pos_to_int)
        .pipe(drop_interaction_csqs)
    )


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
