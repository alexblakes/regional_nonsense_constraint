"""Get phyloP statistics by region."""

import logging
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

import src

FILE_PHYLOP = "data/interim/cds_sites_phylop_pext_missense.tsv"
FILE_REGIONAL_CONSTRAINT = "data/final/regional_nonsense_constraint_for_excel.tsv"
FILE_OUT = "data/interim/phylop_stats_per_region.tsv"
NROWS_PHYLOP = 64_286_505

logger = logging.getLogger(__name__)


def read_phylop(path=FILE_PHYLOP, chunksize=100_000):
    """Read phyloP data in chunks and with a progress bar."""

    n_chunks = (NROWS_PHYLOP / chunksize) + 1

    def reader():
        """Called by tqdm. Separate here for code clarity."""

        return pd.read_csv(
            path,
            sep="\t",
            usecols=["chr", "pos", "enst", "region", "phylop"],
            dtype=defaultdict(lambda: "category", pos="int32", phylop="float32"),
            chunksize=chunksize,
            # nrows=10_000_000,
        )

    # Show a progress bar while reading data
    df = pd.concat(tqdm(reader(), desc="Reading phyloP data.", total=n_chunks))

    logger.info(
        f"phyloP entries: {len(df):,}\n"
        f"Note that phyloP scores are also given for 'transcript' regions."
    )
    logger.info(f"Unique transcripts: {df.enst.nunique():,}")
    logger.info(f"Duplicated by chr, pos: {df.duplicated(['chr','pos']).sum():,}")
    logger.info(
        f"Duplicated by chr, pos, enst, region: "
        f"{df.duplicated(['chr','pos', 'enst', 'region']).sum():,}"
    )

    return df


def get_phylop_stats(df):
    """Get mean, median, and std of phyloP scores per transcript and region."""

    df = (
        df.groupby(["enst", "region"])
        .agg(
            phylop_count=("phylop", "count"),
            phylop_mean=("phylop", "mean"),
            phylop_median=("phylop", "median"),
            phylop_std=("phylop", "std"),
            phylop_sem=("phylop", "sem"),
        )
        .reset_index()
        .astype({"enst": "object", "region": "object"})  # For later merging
    )

    logger.info(f"Entries after aggregation:\n{df.info()}")

    return df


def read_constraint(path=FILE_REGIONAL_CONSTRAINT):
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["symbol", "enst", "region", "oe_ci_hi", "constraint"],
    )

    logger.info(f"Constraint shape: {df.shape}")
    logger.info(f"NaNs in constraint data:\n{df.isna().sum()}")

    return df


def tidy_constraint(df):
    df = df[df["region"] != "transcript"]

    logger.info(f"Constraint shape after dropping 'transcript' regions: {df.shape}")

    return df


def merge_constraint_phylop(left, right):
    df = left.merge(right, how="inner", validate="one_to_one")

    logger.info(f"Shape after merge: {df.shape}")
    logger.info(f"Valid entries after merge:\n{df.notna().sum()}")

    return df


def write_out(df, path=FILE_OUT):
    logger.info(f"Writing to {FILE_OUT}")

    df.to_csv(path, sep="\t", index=False)

    return df


def main():
    """Run as script."""

    phylop = read_phylop().pipe(get_phylop_stats)
    constraint = read_constraint().pipe(tidy_constraint)

    return merge_constraint_phylop(constraint, phylop).pipe(write_out)


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
