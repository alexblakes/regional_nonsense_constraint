"""Get constrained gene lists."""

import logging

import pandas as pd

import src
from src.constraint import get_fdrs_and_gnomad_constraint

_REGIONAL_NONSENSE_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
_GENE_IDS = "data/interim/gene_ids.tsv"
_GNOMAD_V4_CONSTRAINT = "data/raw/gnomad.v4.0.constraint_metrics.tsv"
_FILE_OUT_CANONICAL = "data/final/transcript_list_all.txt"

logger = logging.getLogger(__name__)


def get_transcript_ids(path=_GENE_IDS):
    return pd.read_csv(path, sep="\t").loc[:, "transcript_id"].rename("enst")


def get_constrained_gnomad_transcripts(df):
    m1 = df.pli > 0.9
    m2 = df.loeuf < 0.6

    constrained_transcripts = df[m1 | m2]["enst"].drop_duplicates()

    return constrained_transcripts


def read_regional_nonsense_constraint(path=_REGIONAL_NONSENSE_CONSTRAINT):
    return pd.read_csv(
        path,
        sep="\t",
        usecols=["enst", "region", "constraint"],
    )


def get_constrained_regions(df, region):
    canonical_transcripts = get_transcript_ids()

    constrained_transcripts = (
        df.loc[lambda x: x["region"].isin(region)]
        .query('constraint == "constrained"')
        .loc[lambda x: x.enst.isin(canonical_transcripts), "enst"]
        .drop_duplicates()
    )

    logger.info(
        f'Constrained transcripts in "{region}" region: {len(constrained_transcripts)}'
    )

    return constrained_transcripts


def write_out(series, path):
    series.to_csv(path, sep="\t", index=False)
    return series


def main():
    """Run as script."""

    canonical_transcripts = get_transcript_ids()

    gnomad = (
        get_fdrs_and_gnomad_constraint.get_gnomad_constraint(_GNOMAD_V4_CONSTRAINT)
        .pipe(get_constrained_gnomad_transcripts)
        .loc[lambda x: x.isin(canonical_transcripts)]
    )
    logger.info(f"Constrained transcripts in gnomAD: {len(gnomad)}")

    regional_constraint = read_regional_nonsense_constraint()
    
    any_region = get_constrained_regions(
        regional_constraint, ["nmd_target", "start_proximal", "long_exon", "distal_nmd"]
    )
    nmd_target = get_constrained_regions(regional_constraint, ["nmd_target"])
    start_proximal = get_constrained_regions(regional_constraint, ["start_proximal"])
    long_exon = get_constrained_regions(regional_constraint, ["long_exon"])
    distal = get_constrained_regions(regional_constraint, ["distal_nmd"])

    write_out(canonical_transcripts, _FILE_OUT_CANONICAL)

    transcripts_dict = {
        "gnomad": gnomad,
        "any_region": any_region,
        "nmd_target": nmd_target,
        "start_proximal": start_proximal,
        "long_exon": long_exon,
        "distal": distal,
    }

    for name, transcripts in transcripts_dict.items():
        path = f"data/final/transcript_list_{name}_constrained.txt"
        write_out(transcripts, path)

    return None


if __name__ == "__main__":
    logger = src.setup_logger(__file__)
    main()
