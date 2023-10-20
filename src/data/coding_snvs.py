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


# Functions
def read_getfasta_output(path):
    """Read the output of a bedtools getfasta command from a TSV file."""

    logger.info("Reading getfasta output.")

    df = pd.read_csv(path, sep="\t", header=None, names=["id", "seq"])

    # TODO lose this in production
    _N = 10000
    logger.warning(f"Using only {_N} rows of data for testing.")
    df = df.tail(_N)

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
    assert len(interval) == len(
        sequence
    ), "Interval and sequence are different lengths."

    return list(zip(interval, sequence))


def get_all_sites(df):
    """Get the position and reference allele for sites in all bed intervals."""

    # Get position and reference allele for each site
    logger.info("Getting position and reference allele for each site in the CDS.")
    sites = df.apply(get_pos_and_ref, axis=1)
    sites = sites.explode()
    sites = pd.DataFrame(sites.to_list(), index=sites.index, columns=["pos", "ref"])

    # Manage data types for memory efficiency
    sites["pos"] = pd.to_numeric(sites["pos"], downcast="integer")
    sites = sites.astype(dtype={"ref": "category"})

    logger.info(f"CDS sites: {len(sites)}")

    return sites


def merge_sites_and_bed_intervals(df, sites):
    """Merge sites with bed intervals, on index."""

    # Use categorical dtypes for efficiency
    df["chr"] = df["chr"].astype("category")  
    df = pd.concat([df["chr"], sites], axis=1)
    logger.info(f"CDS sites after merge with bed intervals: {len(df)}")

    # Drop duplicated sites
    df = df.drop_duplicates().sort_values(["chr", "pos"])
    logger.info(f"CDS sites after dropping duplicates: {len(df)}")

    # Sanity checks
    assert (
        df.duplicated(["chr", "pos"]).sum() == 0
    ), "There are sites with multiple reference alleles."
    assert df.isna().sum().sum() == 0, "There are missing values."

    return df


def get_tri_contexts(df):
    """Get the trinucleotide context around each position."""

    # Order positions in each CDS
    df.index = pd.CategoricalIndex(
        df.index, name="cds_id"
    )  # The index acts as a unique ID for each CDS.
    df = df.sort_values(["cds_id", "pos"])

    # Get triplet context
    first = df.groupby(["cds_id"])["ref"].shift(1)  # Order preserved by groupby
    last = df.groupby(["cds_id"])["ref"].shift(-1)
    tri = first.str.cat(others=[df["ref"], last]).rename("tri").astype("category")

    # Merge positions with their trinucleotide contexts
    df = pd.concat([df, tri], axis=1)  # Extreme ends contain NaNs: useful for dropping
    logger.info(
        f"CDS sites before trimming the extended 5' and 3' positions: {len(df)}"
    )

    # CDS intervals were extended by 1 nt each side with the bedtools slop command.
    # This allows annotation of the trinucleotide context around the most 5' and 3'
    # sites. Conveniently, the "tri" column now contains NaNs at the extended positions.
    # These positions are dropped.
    df = df.dropna()
    logger.info(f"CDS sites after trimming: {len(df)}")

    return df


def get_all_possible_alt_alleles(df):
    """Get all possible alt alleles for each position."""

    df = pd.concat([df.copy().assign(alt=base) for base in ["A", "T", "C", "G"]])

    # Convert "alt" column to categorical
    df["alt"] = df["alt"].astype("category")

    # Some reference alleles are "N", so add "N" as a category for alt alleles too
    if "N" in df["ref"].cat.categories:
        df["alt"] = df["alt"].cat.add_categories("N")

    # Remove sites where the reference allele is the same as the alt allele
    df = df[df["ref"] != df["alt"]]

    logger.info(f"Possible SNVs before tidying: {len(df)}")

    return df


def tidy_data(df):
    # Tidy the dataframe

    # Remove sites where the reference allele is "N" (mainly chrY positions)
    df = df[df["ref"] != "N"]
    logger.info(f"Possible SNVs after tidying: {len(df)}")

    return df


def main():
    df = read_getfasta_output(C.CANONICAL_CDS_FASTA).pipe(get_chr_start_end)
    sites = get_all_sites(df)
    df = (
        merge_sites_and_bed_intervals(df, sites)
        .pipe(get_tri_contexts)
        .pipe(get_all_possible_alt_alleles)
        .pipe(tidy_data)
    )

    logger.warning("Lose this return statement when done testing.")
    return df


if __name__ == "__main__":
    main()
    # TODO deal with Ns in the reference sequence

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
