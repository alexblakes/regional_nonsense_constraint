"""Create regional constraint TSV ready for export to Excel."""

import logging

import pandas as pd

import src

FILE_REGIONAL_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
FILE_GENE_IDS = "data/interim/gene_ids.tsv"

logger = logging.getLogger(__name__)


def read_constraint(path=FILE_REGIONAL_CONSTRAINT):
    df = pd.read_csv(path, sep="\t")

    logger.info(f"Constraint entries: {len(df)}")
    logger.info(f"Unique ENSTs in constraint: {df.enst.nunique()}")

    return df


def read_gene_ids(path=FILE_GENE_IDS):
    df = pd.read_csv(path, sep="\t", names=["ensg","enst","symbol"])

    logger.info(f"Gene IDs shape: {df.shape}")
    logger.info(f"Unique ENSGs: {df.ensg.nunique()}")
    logger.info(f"Unique symbols: {df.symbol.nunique()}")

    return df


def merge_constraint_with_ids(left, right):
    df = left.merge(right, how="right", validate="many_to_one")
    logger.info(f"Merged data shape: {df.shape}")
    return df


def main():
    """Run as script."""

    constraint=read_constraint()
    gene_ids=read_gene_ids()

    return merge_constraint_with_ids(constraint, gene_ids)


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
