"""Create regional constraint TSV ready for export to Excel."""

import logging

import pandas as pd

import src

FILE_REGIONAL_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
FILE_GENE_IDS = "data/interim/gene_ids.tsv"
FILE_SHET = "data/raw/genebayes_shet_estimates.tsv"
FILE_OMIM = "data/interim/genemap2_simple.tsv"
FILE_OUT = "data/final/regional_nonsense_constraint_for_excel.tsv"

logger = logging.getLogger(__name__)


def read_constraint(path=FILE_REGIONAL_CONSTRAINT):
    df = pd.read_csv(path, sep="\t")

    logger.info(f"Constraint entries: {len(df)}")
    logger.info(f"Unique ENSTs in constraint: {df.enst.nunique()}")

    return df


def read_gene_ids(path=FILE_GENE_IDS):
    df = pd.read_csv(path, sep="\t", names=["ensg", "enst", "symbol"])

    logger.info(f"Gene IDs shape: {df.shape}")
    logger.info(f"Unique ENSGs: {df.ensg.nunique()}")
    logger.info(f"Unique symbols: {df.symbol.nunique()}")

    return df


def read_shet(path=FILE_SHET):
    df = pd.read_csv(path, sep="\t", usecols=["ensg", "post_mean"]).set_axis(
        ["ensg", "shet_post_mean"], axis=1
    )
    logger.info(f"shet shape: {df.shape}")
    logger.info(f"Unique ENSGs in shet: {df.ensg.nunique()}")
    return df


def read_omim(path=FILE_OMIM):
    df = pd.read_csv(path, sep="\t")
    logger.info(f"OMIM data shape: {df.shape}")
    logger.info(f"Unique ENSGs in OMIM data: {df.ensg.nunique()}")
    return df


def implode_omim(df):
    df = (
        df.groupby("ensg")
        .agg(
            omim_phenotype=("phenotype", ";".join),
            omim_inheritance=("inheritance", ";".join),
        )
        .reset_index()
    )

    logger.info(f"Shape: {df.shape}")
    logger.info(f"Unique ENSGs: {df.ensg.nunique()}")

    return df


def merge_gene_ids(left, right):
    df = left.merge(right, how="inner", validate="many_to_one")
    logger.info(f"Merged data shape: {df.shape}")
    logger.info(f"Unique ENSTs after merge: {df.enst.nunique()}")
    logger.info(
        f"Duplicates by ENSG and region: {df.duplicated(['enst','region']).sum()}"
    )
    logger.info(
        f"Symbols with multiple ENSGs: "
        f"{df[df.duplicated(['symbol','region'])].symbol.unique()}"
    )
    return df


def merge_shet(left, right):
    df = left.merge(right, how="left", validate="many_to_one")
    logger.info(f"Merged data shape: {df.shape}")
    logger.info(f"Unique ENSGs after merge: {df.ensg.nunique()}")
    logger.info(
        f"ENSGs with shet scores: "
        f"{df[['ensg','shet_post_mean']].drop_duplicates('ensg').shet_post_mean.count()}"
    )
    return df


def merge_omim(left, right):
    df = left.merge(right, how="left", validate="many_to_one")

    logger.info(f"Merged data shape: {df.shape}")
    logger.info(f"Unique ENSGs after merge: {df.ensg.nunique()}")
    logger.info(f"ENSGs with OMIM annotations: {df[['ensg','omim_inheritance']].drop_duplicates().omim_inheritance.count()}")
    
    return df


def reorder_columns(df):
    return df[
        [
            "symbol",
            "ensg",
            "enst",
            "region",
            "csq",
            "n_pos",
            "n_obs",
            "n_exp",
            "prop_obs",
            "prop_exp",
            "oe",
            "oe_ci_hi",
            "p",
            "fdr_p",
            "syn_p",
            "constraint",
            "pli",
            "loeuf",
            "gnomad_flags",
            "shet_post_mean",
            "omim_phenotype",
            "omim_inheritance",
        ]
    ]


def write_out(df, path=FILE_OUT):
    logger.info(f"Writing to {FILE_OUT}")
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""

    constraint = read_constraint()
    gene_ids = read_gene_ids()
    shet = read_shet()
    omim = read_omim().pipe(implode_omim)

    df = (
        merge_gene_ids(constraint, gene_ids)
        .pipe(merge_shet, shet)
        .pipe(merge_omim, omim)
        .pipe(reorder_columns)
        .pipe(write_out)
    )

    return df


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
