""" This script finds every possible SNV within the regions of interest.
It also annotates the trinucleotide context around each SNV.

It takes the output of a bedtools getfasta command (tsv) as input.

NB the coordinates of each feature should be extended by 1 each side with
bedtools slop, so that the 3nt context around the most 5' and 3' positions
is still accessible.
"""

# Import relevant modules
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C

# Logging
logger = setup_logger(Path(__file__).stem)

#! TEMP
fasta = C.CANONICAL_CDS_FASTA


# Functions
def read_getfasta_output(path):
    """Read the output of a bedtools getfasta command from a TSV file."""

    logger.info("Reading getfasta output.")

    df = pd.read_csv(path, sep="\t", header=None, names=["id", "seq"])
    logger.warning("Using only five rows of data for testing.")
    df = df.head(5)
    
    logger.info(f"Rows in getfasta output: {len(df)}")

    return df


def get_chr_start_end(df):
    """Extract chromosome, start, and end information from the DataFrame."""

    df[["chr", "coords"]] = df["id"].str.split(":", n=1, expand=True)
    df[["start", "end"]] = df["coords"].str.split("-", n=1, expand=True).astype(int)
    df["start"] = df["start"] + 1  # Revert to 1-based

    df = df[["chr", "start", "end", "seq"]].sort_values(["chr", "start"])

    return df


def get_pos_and_ref(row):
    """Get the position and reference allele for each site in a bed interval."""

    interval = range(row["start"], row["end"] + 1)
    sequence = row["seq"]

    return list(zip(interval, sequence))


def get_all_sites(df):
    """Get the position and reference allele for sites in all bed intervals."""

    logger.info("Getting position and reference allele for each site in the CDS.")
    sites = df.apply(get_pos_and_ref, axis=1)
    sites = sites.explode()

    # Convert sites to a DataFrame
    # Use categorical dtypes for memory efficiency
    sites = pd.DataFrame(
        sites.to_list(), index=sites.index, columns=["pos", "ref"]
    ).astype(dtype={"ref": "category"})

    logger.info(f"CDS sites: {len(sites)}")

    return sites


def merge_sites_and_bed_intervals(df, sites):
    """Merge sites with bed intervals, on index."""

    df["chr"] = df["chr"].astype("category")  # Categorical data for memory efficiency
    df = pd.concat([df["chr"], sites], axis=1)
    logger.info(f"CDS sites after merge with bed intervals: {len(df)}")

    df = df.drop_duplicates()
    logger.info(f"CDS sites after dropping duplicates: {len(df)}")

    assert df.duplicated(["chr","pos"]).sum() == 0, "There are sites with multiple reference alleles."
    assert df.isna().sum().sum() == 0, "There are missing values."

    return df

#! Continue from here


def get_tri_contexts(getfasta_output):
    # Get trinucleotide context around each position.
    df.index = df.index.rename("cds_id")
    df = df.sort_values(["cds_id", "pos"])  # Ensure order of positions

    first = df.groupby(["cds"])["ref"].shift(1)  # Order preserved by groupby
    last = df.groupby(["cds"])["ref"].shift(-1)
    tri = first.str.cat(others=[df["ref"], last]).rename("tri")  # Get triplet context

    df = pd.concat([df, tri], axis=1)  # Extreme ends contain NaNs: useful for dropping
    print(f"They span {len(df)} nt before trimming of the most 3' and 5' positions")

    print(
        f"There are {df.tri.isna().sum()} positions at the extreme ends of the features."
    )
    df = df.dropna()  # Drop extreme end positions
    df = df.sort_values(["chr", "pos"])  # Sort for faster VEP annotation
    print(f"There are {len(df)} nt after trimming.")

    # Get possible alt alleles for each position.
    df["alt"] = [["A", "T", "C", "G"]] * len(df)
    df = df.explode("alt")
    df = df[df["ref"] != df["alt"]]
    df = df[["chr", "pos", "ref", "alt", "tri"]].reset_index(drop=True)
    print(f"There are {df.duplicated().sum()} identical duplicate positions")

    # Tidy the dataframe
    df = df.drop_duplicates()
    df = df[df["ref"] != "N"]  # Mainly chrY positions
    print(f"There are {len(df)} possible SNVs")

    return df


def main():
    df = read_getfasta_output(fasta).pipe(get_chr_start_end)
    sites = get_all_sites(df)
    df = merge_sites_and_bed_intervals(df, sites)
    
    logger.warning("Lose this return statement when done testing.")
    return df


if __name__ == "__main__":
    main()

    # get_fasta_output = "../outputs/gencode_v39_canonical_cds_seq.tsv"

    # # Combine NMD-pos and NMD-esc regions
    # df = get_tri_contexts(get_fasta_output)
    # print(f"There are {len(df)} distinct possible SNVs in all CDS features.")
    # print("Writing trinucleotide contexts to .tsv")
    # df.to_csv("../outputs/cds_trinucleotide_contexts.tsv", sep="\t", index=False)

    # # Create a VCF file for VEP annotation
    # vcf = df.copy().assign(ID=".", QUAL=".", FILTER=".", INFO=".")
    # vcf = vcf[["chr", "pos", "ID", "ref", "alt", "QUAL", "FILTER", "INFO"]]
    # print("Writing SNVs to .vcf")
    # vcf.to_csv(
    #     "../outputs/cds_all_possible_snvs.vcf", sep="\t", index=False, header=False
    # )
