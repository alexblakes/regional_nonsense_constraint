"""Count the number of exons and CDS in each transcript.

Also annotate exon start and end positions. The output is an interim TSV which will be
used to produce NMD annotations for each CDS position.
"""



import pandas as pd
import gtfparse  # * read_gtf makes a call to logging.basicConfig() which overwrites my logging config.

import src
from src import constants as C
from src.regions import canonical_cds as ccds



logger = src.logger


def count_exons(df):
    """Count the number of exons per transcript."""

    df["exon_number"] = df["exon_number"].astype(int)
    df["exon_count"] = df.groupby("transcript_id")["exon_number"].transform("max")

    return df


def get_exon_coords(df):
    """Find exon start and end positions (distinct from CDS starts and ends)."""

    exons = df[df.feature == "exon"]
    exon_coords = exons[["exon_id", "start", "end"]].copy()
    exon_coords.columns = ["exon_id", "exon_start", "exon_end"]
    exon_coords = exon_coords.drop_duplicates()

    logger.info(f"Distinct exons: {len(exon_coords)}")

    return exon_coords


def merge_cds_with_exon_coords(cds, exon_coords):
    """Merge CDS with exon start and end annotations."""

    df = cds.merge(exon_coords, how="left")
    df = df.drop_duplicates(["transcript_id", "exon_id", "start", "end"])

    logger.info(f"Number of CDS after merging with exon coordinates: {len(df)}")
    logger.info(f"Number of CDS positions: {((df.end - df.start) + 1).sum()}")

    return df


def tidy_cds(df):
    """Tidy the CDS dataframe."""

    df = df[
        [
            "seqname",
            "start",
            "end",
            "exon_start",
            "exon_end",
            "strand",
            "transcript_id",
            "exon_id",
            "exon_number",
            "exon_count",
            "cds_number",
            "cds_count",
        ]
    ].copy()

    logger.info(f"Unique transcript IDs: {df['transcript_id'].nunique()}")

    df["transcript_id"] = df["transcript_id"].str.split(".").str[0]

    logger.info(
        f"Unique transcript IDs after dropping version numbers: {df['transcript_id'].nunique()}"
    )
    logger.info(
        f"Duplicated CDS (by exon_id): {df.duplicated(['seqname','exon_id']).sum()}"
    )
    logger.info(
        f"Duplicated CDS (by transcript_id and coords): {df.duplicated(['seqname','transcript_id','start','end']).sum()}"
    )
    logger.info(
        f"Duplicated CDS (by CDS coords): {df.duplicated(['seqname','start','end','strand']).sum()}"
    )
    logger.info(
        f"Duplicated CDS (by exon coords): {df.duplicated(['seqname','exon_start','exon_end','strand']).sum()}"
    )

    return df


def main():
    """Run the script."""

    # Get exons and CDS from GENCODE
    df = (
        gtfparse.read_gtf(C.GENCODE_GTF)
        .pipe(ccds.get_canonical, features=["CDS", "exon"])
        .pipe(count_exons)
    )

    # Annotate CDS with exon counts and coordinates
    exon_coords = get_exon_coords(df)
    cds = df[df["feature"] == "CDS"].copy().pipe(ccds.annotate_cds_number)
    cds = merge_cds_with_exon_coords(cds, exon_coords).pipe(tidy_cds)

    # Write to output
    cds.to_csv(C.CDS_COUNTS_AND_COORDS, sep="\t", index=False)


if __name__ == "__main__":
    src.add_log_handlers()
    main()