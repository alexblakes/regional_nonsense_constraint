"""Get constrained gene lists."""

import logging
from pathlib import Path

import pandas as pd

import src
from src.constraint import constraint_statistics

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_REGIONAL_NONSENSE_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
_GENE_IDS = "data/interim/gene_ids.tsv"
_GNOMAD_V4_CONSTRAINT = "data/raw/gnomad.v4.0.constraint_metrics.tsv"
_FILE_OUT_ALL = "data/final/gene_list_all.txt"
_FILE_OUT_GNOMAD_CST = "data/final/gene_list_gnomad_constrained.txt"
_FILE_OUT_NMD_TARGET = "data/final/gene_list_nmd_target_constrained.txt"
_FILE_OUT_START_PROX = "data/final/gene_list_start_proximal_constrained.txt"
_FILE_OUT_LONG_EXON = "data/final/gene_list_long_exon_constrained.txt"
_FILE_OUT_DISTAL = "data/final/gene_list_distal_constrained.txt"

logger = logging.getLogger(__name__)


def read_gene_ids(path=_GENE_IDS):
    return pd.read_csv(path, sep="\t", usecols=["gene_id", "transcript_id"]).set_axis(
        ["ensg", "enst"], axis="columns"
    )


def get_constrained_gnomad_transcripts(gnomad):
    logger.info("Getting constrained gnomAD transcripts.")

    m1 = gnomad.pli > 0.9
    m2 = gnomad.loeuf < 0.6

    gnomad_strong = gnomad[m1 | m2]["enst"].drop_duplicates()

    logger.info(f"Constrained transcripts in gnomAD: {len(gnomad_strong)}")

    return gnomad_strong


def get_gene_ids(df, gene_ids):
    """Get gene IDs for a dataframe with a list of transcript IDs."""

    df = gene_ids.merge(df, how="inner", on="enst", validate="one_to_one")["ensg"]
    assert df.duplicated().sum() == 0, "Duplicated gene IDs."

    logger.info(f"Gene ids: {len(df)}")

    return df


def get_regional_constrained_genes(
    path=_REGIONAL_NONSENSE_CONSTRAINT, region="nmd_target"
):
    logger.info(region)

    rnc = (
        pd.read_csv(
            path,
            sep="\t",
            usecols=["enst", "region", "constraint"],
        )
        .query("region == @region")
        .query('constraint == "constrained"')
    )

    assert (
        rnc.enst.duplicated().sum() == 0
    ), f"Duplicated transcript IDs in {region} region."

    assert len(rnc) > 0, "Empty dataframe."

    return rnc


def main():
    """Run as script."""

    # Read data
    gene_ids = read_gene_ids(_GENE_IDS)
    gnomad = constraint_statistics.get_gnomad_constraint(_GNOMAD_V4_CONSTRAINT)

    # Define gene lists
    ## All genes with a protein-coding canonical transcript
    _all = gene_ids["ensg"].drop_duplicates()
    logger.info(f"All genes: {len(_all)}")

    ## Constrained genes in gnomAD
    gnomad_strong = get_constrained_gnomad_transcripts(gnomad).pipe(
        get_gene_ids, gene_ids
    )

    ## Genes with regional nonsense constraint
    nmd_target = get_regional_constrained_genes().pipe(get_gene_ids, gene_ids)
    start_proximal = get_regional_constrained_genes(region="start_proximal").pipe(
        get_gene_ids, gene_ids
    )
    long_exon = get_regional_constrained_genes(region="long_exon").pipe(
        get_gene_ids, gene_ids
    )
    distal = get_regional_constrained_genes(region="distal_nmd").pipe(
        get_gene_ids, gene_ids
    )

    # Write to output
    write_out = lambda x, y: x.to_csv(y, index=False, header=None)

    write_out(_all, _FILE_OUT_ALL)
    write_out(gnomad_strong, _FILE_OUT_GNOMAD_CST)
    write_out(nmd_target, _FILE_OUT_NMD_TARGET)
    write_out(start_proximal, _FILE_OUT_START_PROX)
    write_out(long_exon, _FILE_OUT_LONG_EXON)
    write_out(distal, _FILE_OUT_DISTAL)

    return None


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
