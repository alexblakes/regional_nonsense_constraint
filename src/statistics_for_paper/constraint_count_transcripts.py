"""Count transcripts with regional nonsense constraint."""



import pandas as pd

import src
from src.constraint import get_fdrs_and_gnomad_constraint as fgc
from src.constraint import gene_lists as gl

_REGIONAL_NONSENSE_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
_GNOMAD_V4_CONSTRAINT = "data/raw/gnomad.v4.0.constraint_metrics.tsv"
_GENE_IDS = "data/interim/gene_ids.tsv"
logger = src.logger


def main():
    """Run as script."""
    enst_ids = gl.read_gene_ids().loc[:, "enst"].pipe(set)

    gnomad = (
        fgc.get_gnomad_constraint(_GNOMAD_V4_CONSTRAINT)
        .query("pli > 0.9 | loeuf < 0.6")
        .loc[:, "enst"]
        .pipe(set)
    )

    nmd_target = (
        pd.read_csv(_REGIONAL_NONSENSE_CONSTRAINT, sep="\t")
        .query("region == 'nmd_target'")
        .query("constraint == 'constrained'")
        .loc[:, "enst"]
        .pipe(set)
        .intersection(enst_ids)
    )

    logger.info(
        f"Constrained transcripts in gnomAD: {len(gnomad.intersection(enst_ids))}"
    )

    return nmd_target


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
