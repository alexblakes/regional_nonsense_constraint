"""Module docstring"""

# # Canonical CDS
# This script defines the coordinates of the CDS for canonical transcript in a
# GENCODE annotation.
# It takes a GENCODE .gtf file as input.
# It gives a .bed output for each of the desired features.

# Import modules
import numpy as np
import pandas as pd
import gtfparse

dx download -f -o ../data/ data/gencode.v39.annotation.gtf


def get_gencode_gtf(path):
    """Read a GENCODE .gtf into memory with gtfparse"""
    gtf = gtfparse.read_gtf(path)
    return gtf

def get_canonical_cds(gtf):
    """Identify all CDS features in each Ensembl_canonical from GENCODE"""

    # Subset to Ensembl_canonical CDS features in protein coding genes
    canonical_cds = gtf[
        (gtf.feature == "CDS")
        & (gtf.tag.str.contains("Ensembl_canonical"))
        & (gtf.gene_type == "protein_coding")
    ].copy()

    return canonical_cds


def annotate_exon_number(cds):
    """Count the number of CDS exons in each transcript"""
    cds["exon_number"] = cds["exon_number"].astype(int)
    cds["cds_number"] = cds.groupby("transcript_id")["exon_number"].rank().astype(int)
    return cds
def gtf_to_bed(gtf, id_list):
    """Reformat a .gtf to .bed format.
    Any identifiers (.gtf file column names) given in "id_list", will be included as
    a comma-separated string in the bed id column.
    """
    bed = gtf.copy()

    # Give desired identifiers in the .bed id column, as a comma-separated string
    bed["id"] = bed[id_list].apply(
        lambda row: ",".join(row.values.astype("str")),
        axis=1,
    )

    # Reformat to .bed
    bed.loc[:, "score"] = "."
    bed.loc[:, "start"] = bed["start"] - 1
    bed = (
        bed[["seqname", "start", "end", "id", "score", "strand"]]
        .copy()
        .sort_values(by=["seqname", "start"])
        .drop_duplicates()
    )

    return bed

def write_bed(bed, path, chr_prefix="chr"):
    """Write a .bed file to output.
    Specify the chr_prefix; either "chr" (GRCh38) or "" (GRCh37).
    """
    # Add or remove the "chr" prefix
    bed = bed.copy()
    bed["seqname"] = chr_prefix + bed["seqname"].str.slice(start=3)

    # Write to output
    bed.to_csv(path, sep="\t", header=False, index=False)

    return bed


if __name__ == "__main__":
    # Read GTF data
    gencode_path = "../data/gencode.v39.annotation.gtf"
    gtf = get_gencode_gtf(gencode_path)

    # Define regions of interest
    cds = get_canonical_cds(gtf).pipe(annotate_exon_number)

    # Save intermediate outputs
    ## Canonical transcript list
    pd.Series(cds.transcript_id.unique()).to_csv(
        "../outputs/canonical_transcripts.txt", sep="\t", index=False, header=False
    )
    ## ENSG, ENST, and HGNC symbols
    cds[["gene_id", "transcript_id", "gene_name"]].drop_duplicates().to_csv(
        "../outputs/gene_ids.tsv", sep="\t", index=False
    )

    # Convert .gtf data to .bed format
    bed_ids = ["gene_id", "transcript_id", "exon_id", "cds_number"]
    cds_bed = gtf_to_bed(cds, bed_ids)

    # Write to output
    gencode_version = "v39"
    feature = "cds"
    out_path_chr = f"../outputs/gencode_{gencode_version}_canonical_{feature}_chr.bed"

    write_bed(cds_bed, out_path_chr, "chr")


