"""
Get CDS features from canonical transcripts in a GENCODE GTF.

This module extracts the coding sequence (CDS) for each canonical transcript in the 
input GTF file, annotates annotates the number of CDS exons in each transcript, and 
writes the results to a BED file.

Functions:
    - get_gencode_gtf(path): Reads a GENCODE .gtf into memory with gtfparse.
    - get_canonical_cds(gtf): Subsets to Ensembl_canonical CDS features in protein 
        coding genes.
    - annotate_exon_number(gtf): Counts the number of CDS exons in each transcript.
    - gtf_to_bed(gtf, ids): Converts a .gtf file to .bed format.
    - write_bed(bed, path, chr_prefix="chr"): Writes a .bed file to output.
    - main(): Runs all functions in the moduule.
"""


# Imports
from pathlib import Path

import pandas as pd
import gtfparse # * read_gtf makes a call to logging.basicConfig() which overwrites my logging config.

from src import constants as C
from src import setup_logger


# Module constants
_BED_COLUMNS = ["seqname", "start", "end", "id", "score", "strand"]
_BED_IDS = ["gene_id", "transcript_id", "exon_id", "cds_number"]
_CHR_PREFIX = "chr"


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_canonical(gtf, features=["CDS"]):
    """Subset to features in Ensembl_canonical protein coding genes."""

    logger.info("Filtering for canonical features...")

    cds = gtf[
        (gtf.feature.isin(features))
        & (gtf.tag.str.contains("Ensembl_canonical"))
        & (gtf.gene_type == "protein_coding")
    ].copy()

    logger.info(f"Chromosomes represented: {cds.seqname.unique()}")
    logger.info(f"Unique gene IDs: {cds.gene_id.nunique()}")
    logger.info(f"Unique transcript IDs: {cds.transcript_id.nunique()}")
    logger.info(f"Feature counts:\n{cds.feature.value_counts()}")

    return cds


def annotate_cds_number(gtf):
    """Count the number of CDS exons in each transcript."""

    gtf["exon_number"] = gtf["exon_number"].astype(int)
    gtf["cds_number"] = gtf.groupby("transcript_id")["exon_number"].rank(method="min").astype(int)
    gtf["cds_count"] = gtf.groupby("transcript_id")["cds_number"].transform("max")

    return gtf


def gtf_to_bed(gtf, ids):
    """
    Convert a .gtf file to .bed format.

    Args:
        gtf (pandas.DataFrame): A pandas DataFrame in gtf format.
        ids (list): A list of identifiers to include in the .bed id column.

    Returns:
        pandas.DataFrame: A pandas DataFrame in .bed format.
    """

    logger.info("Converting to bed...")

    # Place a comma-separated string of identifiers in the .bed id column
    gtf["id"] = gtf[ids].apply(
        lambda row: ",".join(row.values.astype("str")),
        axis=1,
    )

    # Reformat to .bed
    gtf.loc[:, "score"] = "."
    gtf.loc[:, "start"] = gtf["start"] - 1
    bed = (
        gtf[_BED_COLUMNS].copy().sort_values(by=["seqname", "start"]).drop_duplicates()
    )

    logger.info(f"Number of unique IDs in the bed file: {bed.id.nunique()}")
    logger.info(
        f"Number of unique bed intervals: {bed.duplicated(['seqname','start','end']).value_counts()[False]}"
    )

    return bed


def write_bed(bed, path, chr_prefix="chr"):
    """
    Write a .bed file to output.

    Args:
    - bed (pandas.DataFrame): A DataFrame in bed format.
    - path (str): Output file path.
    - chr_prefix (str): The prefix to add to the seqname column. Default is "chr".

    Returns:
    - None
    """

    logger.info(f"Writing bed to {path}")

    # Add or remove the "chr" prefix
    bed["seqname"] = chr_prefix + bed["seqname"].str.slice(start=3)

    # Write to output
    bed.to_csv(path, sep="\t", header=False, index=False)

    return None


def main():
    """Runs all functions in this module."""

    (
        gtfparse.read_gtf(C.GENCODE_GTF, features="CDS")
        .pipe(get_canonical)
        .pipe(annotate_cds_number)
        .pipe(gtf_to_bed, _BED_IDS)
        .pipe(write_bed, C.CANONICAL_CDS_BED, _CHR_PREFIX)
    )

    return None


if __name__ == "__main__":
    main()
