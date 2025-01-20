"""Extract Zoonomia phyloP scores for 241 mammalian species."""

# Imports

import pandas as pd
import pyBigWig

from src import constants as C
from src import setup_logger


# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def read_cds_bed(path):
    """Read the canonical CDS bedfile to memory."""

    logger.info("Reading CDS bed file.")

    bed = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "id", "score", "strand"],
        usecols=["chr", "start", "end"],
        # nrows=100,  #! Testing
    )

    return bed


def annotate_phylop_scores(bed, bw):
    """Annotate CDS sites with phyloP scores."""

    #! chrM sites are dropped; they are "out of bounds" in the phyloP bw file.
    logger.info("Dropping mitochondrial sites.")
    bed = bed.query("chr != 'chrM'").copy()

    logger.info("Annotating CDS sites with phyloP scores.")

    bed["phylop"] = bed.apply(
        lambda x: bw.intervals(x["chr"], x["start"], x["end"]), axis=1
    )

    bed = bed.drop(["start", "end"], axis=1).explode("phylop").dropna().set_index("chr")

    bed = (
        pd.DataFrame(
            [[*a] for a in bed.phylop.values],
            columns=["start", "pos", "phylop"],
            index=bed.index,
        )
        .reset_index(drop=False)
        .drop("start", axis=1)
        .drop_duplicates()
    )

    logger.info(f"CDS sites with phyloP annotation: {len(bed)}")
    logger.info(f"Duplicated sites: {bed.duplicated(['chr','pos']).sum()}")

    return bed


def main():
    """Run as script."""

    bw = pyBigWig.open(C.PHYLOP_BW)
    bed = read_cds_bed(C.CANONICAL_CDS_BED)

    df = annotate_phylop_scores(bed, bw)

    # Write to output
    logger.info("Writing to output.")
    df.to_csv(C.PHYLOP_CDS_SCORES, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    main()
