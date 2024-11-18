"""Merge regional, shet, and gnomAD constraint annotations."""

import logging

import pandas as pd
from scipy import stats

import src

logger = logging.getLogger(__name__)

FILE_REGIONAL_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
FILE_GENE_IDS = "data/interim/gene_ids.tsv"
FILE_SHET = "data/raw/genebayes_shet_estimates.tsv"
FILE_OUT = "data/interim/shet_gnomad_regional_constraint.tsv"


def read_regional_annotation(path=FILE_REGIONAL_CONSTRAINT):
    df = pd.read_csv(
        path, sep="\t", usecols=["enst", "region", "constraint", "pli", "loeuf"]
    )

    logger.info(f"Regional constraint entries: {len(df)}")
    logger.info(
        f"Constraint by region:\n{df.groupby('region').constraint.value_counts(dropna=False)}"
    )

    return df


def reshape_regional(df):
    df = df.pivot(
        index=["enst", "pli", "loeuf"], columns="region", values="constraint"
    ).reset_index()

    logger.info(f"Entries after reshaping: {len(df)}")
    logger.info(f"Unique ENST IDs: {df.enst.nunique()}")

    return df


def make_regional_boolean(df):
    df = df.replace(
        {"constrained": True, "unconstrained": False, "indeterminate": False}
    )
    regions = ["transcript", "nmd_target", "start_proximal", "long_exon", "distal_nmd"]

    df.loc[:, regions] = df.loc[:, regions].fillna(False)

    return df


def read_gene_ids(path=FILE_GENE_IDS):
    df = pd.read_csv(path, sep="\t", names=["ensg", "enst", "symbol"]).drop(
        "symbol", axis=1
    )

    logger.info(f"Gene IDs length: {len(df)}")
    logger.info(f"Unique ENSGs: {df.ensg.nunique()}")
    logger.info(f"Unique ENSTs: {df.enst.nunique()}")

    return df


def merge_gene_ids(left, right, **kwargs):
    logger.info(f"Length before merge: {len(left)}")

    df = left.merge(right, **kwargs)

    logger.info(f"Length after merge: {len(df)}")
    logger.info(f"Unique ENSGs: {df.ensg.nunique()}")

    return df


def read_shet(path=FILE_SHET):
    df = pd.read_csv(path, sep="\t", usecols=["ensg", "post_mean"]).rename(
        columns={"post_mean": "shet_post_mean"}
    )
    logger.info(f"Length shet data: {len(df)}")
    logger.info(f"Unique ENSGs in shet data: {df.ensg.nunique()}")

    return df


def find_constrained_genes(df):
    quantile = (
        stats.percentileofscore(df.loeuf, 0.6, kind="strict", nan_policy="omit") / 100
    )
    shet_cutoff = df.shet_post_mean.quantile(1 - quantile)

    logger.info(f"Proportion of LOEUF scores < 0.6: {quantile}")
    logger.info(f"shet cutoff at {1 - quantile:2f} quantile: {shet_cutoff}")

    df = (
        df.assign(
            any_region_constrained=lambda x: x.loc[
                :, ["distal_nmd", "long_exon", "nmd_target", "start_proximal"]
            ].any(axis=1)
        )
        .assign(shet_constrained=lambda x: x.shet_post_mean > shet_cutoff)
        .assign(gnomad_constrained=lambda x: (x.loeuf < 0.6) | (x.pli > 0.9))
    )

    logger.info(f"Constrained in any region: {df.any_region_constrained.sum()}")
    logger.info(f"gnomAD constrained: {df.gnomad_constrained.sum()}")
    logger.info(f"shet constrained: {df.shet_constrained.sum()}")

    return df


def reorder_columns(df):
    return df[
        [
            "ensg",
            "transcript",
            "nmd_target",
            "start_proximal",
            "long_exon",
            "distal_nmd",
            "any_region_constrained",
            "gnomad_constrained",
            "shet_constrained",
        ]
    ]


def write_out(df, path=FILE_OUT):
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""

    gene_ids = read_gene_ids()
    shet = read_shet()

    regional = (
        (read_regional_annotation().pipe(reshape_regional).pipe(make_regional_boolean))
        .pipe(merge_gene_ids, gene_ids, how="inner", validate="one_to_one")
        .pipe(merge_gene_ids, shet, how="left", validate="one_to_one")
        .pipe(find_constrained_genes)
        .pipe(reorder_columns)
        .pipe(write_out)
    )

    #! 18647 -> 18369 transcripts after merging with gene ids.

    return regional


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
