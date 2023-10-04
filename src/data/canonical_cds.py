"""Canonical CDS

This script defines the coordinates of the CDS for canonical transcripts in a GENCODE
annotation. It takes a GENCODE .gtf file as input. It writes a .bed output for each of 
the desired features.

Functions:
- get_gencode_gtf(path): Reads a GENCODE .gtf into memory with gtfparse.
- get_canonical_cds(gtf): Subsets to Ensembl_canonical CDS features in protein coding genes.
- annotate_exon_number(cds): Counts the number of CDS exons in each transcript.
- gtf_to_bed(gtf, id_list): Converts a .gtf file to .bed format.
- write_bed(bed, path, chr_prefix="chr"): Writes a .bed file to output.
- canonical_cds(in_file, out_file): Extracts the coding sequence (CDS) for each canonical transcript in the input GTF 
file, annotates each CDS with its exon number, and writes the results to a BED file.
- parse_args(): Parses command line arguments.
- main(): Calls the canonical_cds function with parsed arguments.

Usage:
- python3 canonical_cds.py --in_file <path_to_input_gtf> --out_file <path_to_output_bed>

Returns:
- None
"""

# Import modules
import pandas as pd
import gtfparse
import argparse


# Functions
def get_gencode_gtf(path):
    """Read a GENCODE .gtf into memory with gtfparse."""

    return gtfparse.read_gtf(path)


def get_canonical_cds(gtf):
    """Subset to Ensembl_canonical CDS features in protein coding genes."""

    canonical_cds = gtf[
        (gtf.feature == "CDS")
        & (gtf.tag.str.contains("Ensembl_canonical"))
        & (gtf.gene_type == "protein_coding")
    ]

    return canonical_cds


def annotate_exon_number(cds):
    """Count the number of CDS exons in each transcript."""

    cds["exon_number"] = cds["exon_number"].astype(int)
    cds["cds_number"] = cds.groupby("transcript_id")["exon_number"].rank().astype(int)

    return cds


def gtf_to_bed(gtf, id_list):
    """
    Convert a .gtf file to .bed format.

    Args:
        gtf (pandas.DataFrame): A pandas DataFrame in gtf format.
        id_list (list): A list of identifiers to include in the .bed id column.

    Returns:
        pandas.DataFrame: A pandas DataFrame in .bed format.
    """
    
    # Give desired identifiers in the .bed id column, as a comma-separated string
    gtf["id"] = gtf[id_list].apply(
        lambda row: ",".join(row.values.astype("str")),
        axis=1,
    )

    # Reformat to .bed
    gtf.loc[:, "score"] = "."
    gtf.loc[:, "start"] = gtf["start"] - 1
    bed = (
        gtf[["seqname", "start", "end", "id", "score", "strand"]]
        .copy()
        .sort_values(by=["seqname", "start"])
        .drop_duplicates()
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

    # Add or remove the "chr" prefix
    bed = bed.copy()
    bed["seqname"] = chr_prefix + bed["seqname"].str.slice(start=3)

    # Write to output
    bed.to_csv(path, sep="\t", header=False, index=False)

    return None


def canonical_cds(in_file, out_file):
    """
    Extracts the coding sequence (CDS) for each canonical transcript in the input GTF 
    file, annotates each CDS with its exon number, and writes the results to a BED file.

    Args:
        in_file (str): Path to the input GTF file.
        out_file (str): Path to the output BED file.

    Returns:
        None
    """

    bed_ids = ["gene_id", "transcript_id", "exon_id", "cds_number"]

    (
        get_gencode_gtf(in_file)
        .pipe(get_canonical_cds)
        .pipe(annotate_exon_number)
        .pipe(gtf_to_bed, bed_ids)
        .pipe(write_bed, out_file)
    )

    return None


def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument(
        "--in_file",
        type=str,
        default="data/external/gencode.v39.annotation.gtf",
        help="Path to GENCODE .gtf file",
    )
    parser.add_argument(
        "--out_file",
        type=str,
        default="data/interim/gencode_v39_canonical_cds.bed",
        help="Path to output .bed file",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    canonical_cds(args)

    return None


if __name__ == "__main__":
    main()
