"""
Assign NMD annotations to all CDS positions.

Label all CDS positions with an NMD annotation. The annotations include:
1) Start-proximal (<150nt from translation start site)
2) Long exons (>400nt upstream of the splice donor site)
3) Last exon
4) 50nt rule (within the most 3' 50nt of the penultimate exon)
"""

# Imports
from pathlib import Path

import pandas as pd
import numpy as np

from src import setup_logger
from src import constants as C

# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def get_cds_positions(df):
    """Get positions of every nucleotide in the CDS"""

    df = df.set_index(["seqname", "transcript_id", "exon_id"])
    df["pos"] = df.apply(lambda x: list(range(x["start"], x["end"] + 1)), axis=1)
    df = df.explode("pos").astype({"pos": int})

    return df


def get_feature_lengths(df):
    """Get the length of the CDS and CDS exons."""

    df["cds_len"] = df.groupby(level=["seqname", "transcript_id"])["pos"].transform(
        "count"
    )
    df["cds_exon_len"] = df.groupby(level="exon_id")["pos"].transform("count")

    # Sanity checks
    assert (df["cds_exon_len"] <= 0).sum() == 0
    assert (df["cds_len"] <= 0).sum() == 0

    return df


def annotate_start_distance(df, strand):
    """Annotate the distance of CDS positions from the translation start site."""

    if strand == "+":
        ascending = True
    if strand == "-":
        ascending = False

    df["start_distance"] = (
        df.groupby(level=["seqname", "transcript_id"])["pos"]
        .rank(ascending=ascending)
        .astype(int)
    )

    # Sanity checks
    assert (df["start_distance"] <= 0).sum() == 0

    return df


def annotate_splice_donor_distance(df, strand):
    """Annotate the distance of CDS positions from the downstream splice donor site."""

    if strand == "+":
        _5_prime = df["pos"]
        _3_prime = df["exon_end"]
    if strand == "-":
        _5_prime = df["exon_start"]
        _3_prime = df["pos"]

    # Find the distance of each CDS position to the splice donor site
    df["splice_donor_distance"] = (_3_prime - _5_prime) + 1

    # Where the CDS is in the last exon, there is no downstream splice donor site
    df.loc[df["exon_number"] == df["exon_count"], "splice_donor_distance"] = np.nan

    # Sanity check
    assert (df["splice_donor_distance"] <= 0).sum() == 0

    return df


def annotate_last_exons(df):
    """Annotate CDS positions in a last exon."""

    df["last_exon"] = np.where(df["exon_count"] == df["exon_number"], 1, 0)

    return df


def annotate_fifty_nt_rule(df, strand):
    """Annotate CDS positions subject to the 50nt rule."""

    if strand == "+":
        _5_prime = df["pos"]
        _3_prime = df["exon_end"]
    if strand == "-":
        _5_prime = df["exon_start"]
        _3_prime = df["pos"]

    df["fifty_nt_rule"] = np.where(
        (df["exon_number"] == df["exon_count"] - 1) & ((_3_prime - _5_prime) + 1 <= 50),
        1,
        0,
    )

    return df


def get_nmd_criteria(df, strand):
    """Annotate positions with all NMD criteria."""

    df = (
        annotate_start_distance(df, strand)
        .pipe(annotate_splice_donor_distance, strand)
        .pipe(annotate_last_exons)
        .pipe(annotate_fifty_nt_rule, strand)
    )

    return df


def get_nmd_annotations(df):
    """Get NMD annotations for each site."""

    logger.info(f"Start proximal sites: {(df['start_distance'] < 150).sum()}")
    logger.info(f"Long exon sites: {(df['splice_donor_distance'] > 400).sum()}")
    logger.info(f"Last exon sites: {df['last_exon'].sum()}")
    logger.info(f"Fifty nt rule sites: {df['fifty_nt_rule'].sum()}")

    a = pd.Series(np.where(df["start_distance"] <= 150, "start_proximal,", ""))
    b = pd.Series(np.where(df["splice_donor_distance"] > 400, "long_exon,", ""))
    c = pd.Series(np.where(df["fifty_nt_rule"] == 1, "fifty_nt,", ""))
    d = pd.Series(np.where(df["last_exon"] == 1, "last_exon,", ""))

    # Where sites have overlapping NMD annotations, join them
    df["nmd"] = pd.Series(["".join([w, x, y, z]) for w, x, y, z in zip(a, b, c, d)])

    # Sites with no NMD-escape annotation are NMD targets
    df["nmd"] = df["nmd"].replace("", "nmd_target")

    logger.info(f"Overlapping NMD annotations value counts:\n{df.nmd.value_counts()}")

    return df


def get_definitive_nmd_annotations(df):
    """Give a definitive NMD annotation to each site.

    "last_exon" and "fifty_nt" annotations are merged into "distal_nmd".
    Overlapping annotations are removed, with this priority:
    "start_proximal" > "distal_nmd" > "long_exon"
    """

    df["nmd_definitive"] = df["nmd"].copy()

    # Start-proximal sites
    df.loc[
        df["nmd_definitive"].str.contains("start_proximal"), "nmd_definitive"
    ] = "start_proximal"

    # Distal NMD sites
    df.loc[
        (df["nmd_definitive"].str.contains("fifty_nt"))
        | (df["nmd_definitive"].str.contains("last_exon")),
        "nmd_definitive",
    ] = "distal_nmd"

    # Long exon sites
    df["nmd_definitive"] = df["nmd_definitive"].replace("long_exon,", "long_exon")

    logger.info(f"NMD annotation value counts:\n{df.nmd_definitive.value_counts()}")

    return df


def tidy_dataframe(df):
    """Tidy the data and add miscellaneous annotations."""

    # Get relative CDS position
    df["relative_cds_position"] = df["start_distance"] / df["cds_len"]

    # Select columns
    df = (
        df[
            [
                "seqname",
                "pos",
                "transcript_id",
                "exon_id",
                "strand",
                "start",
                "end",
                "exon_number",
                "exon_count",
                "cds_number",
                "cds_count",
                "cds_len",
                "cds_exon_len",
                "relative_cds_position",
                "start_distance",
                "splice_donor_distance",
                "last_exon",
                "fifty_nt_rule",
                "nmd",
                "nmd_definitive",
            ]
        ]
    ).rename(columns={"seqname": "chr"})

    logger.info(f"Number of CDS positions:{len(df)}")
    logger.info(
        f"Number of unique CDS positions:{len(df.drop_duplicates(['chr','pos']))}"
    )

    # Sanity checks
    assert (df.duplicated(["chr","pos","transcript_id"]).sum()) == 0

    return df


def main():
    """Run the script."""

    # Get positions within CDS
    cds = (
        pd.read_csv(C.CDS_COUNTS_AND_COORDS, sep="\t")
        .pipe(get_cds_positions)
        .pipe(get_feature_lengths)
    )

    # Annotate NMD criteria for positions in fwd and rev trancripts
    fwd = cds[cds["strand"] == "+"].copy().pipe(get_nmd_criteria, strand="+")
    rev = cds[cds["strand"] == "-"].copy().pipe(get_nmd_criteria, strand="-")

    # Combine fwd and rev transcripts, get NMD annotations
    df = (
        pd.concat([fwd, rev])
        .reset_index()
        .pipe(get_nmd_annotations)
        .pipe(get_definitive_nmd_annotations)
        .pipe(tidy_dataframe)
    )

    # Write to output
    df.to_csv(C.NMD_ANNOTATIONS, sep="\t", index=False)

    return df  # TODO Testing


if __name__ == "__main__":
    main()
