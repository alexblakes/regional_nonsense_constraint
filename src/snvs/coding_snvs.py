""" This script finds every possible SNV within the regions of interest.
It also annotates the trinucleotide context around each SNV.

It takes the output of a bedtools getfasta command (tsv) as input.
"""

import logging
from pathlib import Path

import pandas as pd
import numpy as np

import src
from src import constants as C

_CHUNKS = 30
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

logger = logging.getLogger(__name__)


def read_getfasta_output(path):
    """Read the output of a bedtools getfasta command from a TSV file."""

    df = pd.read_csv(path, sep="\t", header=None, names=["id", "seq"])
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
    assert len(interval) == len(sequence)

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

    logger.info(f"CDS sites (raw): {len(sites)}")

    return sites


def merge_sites_and_bed_intervals(df, sites):
    """Merge sites with bed intervals, on index."""

    # Use categorical dtypes for efficiency
    df["chr"] = df["chr"].astype("category")

    # Merge sites with bed intervals
    df = pd.concat([df["chr"], sites], axis=1)

    # Don't drop duplicated sites yet. Annotation of the trinucleotide context depends
    # on the sites within each CDS having monotonically increasing values.

    # Sanity checks
    assert df.isna().sum().sum() == 0, "There are missing values."

    return df


def get_tri_contexts(df):
    """Get the trinucleotide context around each position."""

    # Order positions in each CDS
    df.index = pd.CategoricalIndex(df.index, name="cds_id")
    df = df.sort_values(["cds_id", "pos"])

    # Group by CDS feature
    df_grouped = df.groupby("cds_id")

    # Check that the position in each interval increases by 1 each time
    assert (
        df_grouped["pos"].diff().fillna(1) == 1
    ).all(), "Positions are not adjacent."

    # Get triplet context
    first = df_grouped["ref"].shift(1)  # Order preserved by groupby
    last = df_grouped["ref"].shift(-1)
    tri = first.str.cat(others=[df["ref"], last]).rename("tri").astype("category")

    # Merge positions with their trinucleotide contexts
    df = pd.concat([df, tri], axis=1)  # Extreme ends contain NaNs: useful for dropping

    return df


def tidy_tri_contexts(df):
    """Tidy the dataframe of trinucleotide contexts.

    Drop duplicated sites, sites with "N" as the reference allele, and the extreme 5'
    and 3' positions adjacent to the CDS.

    CDS intervals had been extended by 1 nt in
    both 5' and 3' directions with the bedtools slop command. This allows annotation of
    the trinucleotide context around the most 5' and 3' CDS sites. Conveniently, the
    "tri" column now contains NaNs at the extended positions. These positions are
    dropped.
    """

    # Drop extended 5' and 3' positions
    df = df.dropna()
    logger.info(f"CDS sites after trimming extended 5' and 3' positions: {len(df)}")

    # Drop sites where the reference allele is "N"
    df = df[df["ref"] != "N"]
    df["ref"] = df["ref"].cat.remove_unused_categories()
    logger.info(f"CDS sites after dropping 'N' reference alleles: {len(df)}")

    # Drop duplicated sites
    df = df.drop_duplicates().sort_values(["chr", "pos"])
    assert df.duplicated(["chr", "pos"]).sum() == 0, "There are duplicated sites."
    logger.info(f"CDS sites after dropping duplicates: {len(df)}")

    return df


def get_all_possible_alt_alleles(df):
    """Get all possible alt alleles for each position."""

    # Reset index for later sorting
    df = df.sort_values(["chr", "pos"]).reset_index(drop=True)

    # Get all possible alt alleles
    df = pd.concat([df.copy().assign(alt=base) for base in ["A", "T", "C", "G"]])

    # Convert "alt" column to categorical
    df["alt"] = df["alt"].astype("category")

    # Remove sites where the reference allele is the same as the alt allele
    df = df[df["ref"] != df["alt"]]
    logger.info(f"Possible SNVs: {len(df)}")

    # Sort by chromosome and position for use with VEP
    df = df[["chr", "pos", "ref", "alt", "tri"]].sort_index()

    return df


def convert_to_vcf(df):
    """Convert a DataFrame to VCF format."""

    for col in ["ID", "QUAL", "FILTER", "INFO"]:
        df[col] = pd.Categorical(["."] * len(df))

    df = df[["chr", "pos", "ID", "ref", "alt", "QUAL", "FILTER", "INFO"]]

    return df


def main():
    """Run the script."""

    df = read_getfasta_output(C.CANONICAL_CDS_FASTA).pipe(get_chr_start_end)
    sites = get_all_sites(df)
    df = (
        merge_sites_and_bed_intervals(df, sites)
        .pipe(get_tri_contexts)
        .pipe(tidy_tri_contexts)
        .pipe(get_all_possible_alt_alleles)
    )

    # Save .tsv with all trinucleotide contexts
    logger.info("Writing TSV")
    df.to_csv(C.CDS_ALL_SNVS_TRI_CONTEXT, sep="\t", index=False)

    # Save to VCF for VEP annotation
    logger.info("Writing VCF")
    df = convert_to_vcf(df)
    df.to_csv(C.CDS_ALL_SNVS_VCF, sep="\t", index=False, header=False)

    # Split output for VEP annotation
    logger.info("Splitting output for parallel VEP annotation")
    for idx, chunk in enumerate(np.array_split(df, _CHUNKS)):
        chunk.to_csv(
            f"{C.VEP}/in_{idx}.vcf",
            sep="\t",
            index=False,
            header=False,
        )

    logger.info("Done.")


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()