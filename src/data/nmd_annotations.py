"""
Assign NMD annotations to all CDS positions

Label all CDS positions with an NMD annotation. The annotations include:
1) Start-proximal (<150nt from translation start site)
2) Long exons (>400nt upstream of the splice donor site)
3) Last exon
4) 50nt rule (within the most 3' 50nt of the penultimate exon)
"""

# Imports
from pathlib import Path

import pandas as pd
import gtfparse  # * read_gtf makes a call to logging.basicConfig() which overwrites my logging config.

from src import setup_logger
from src import constants as C
from src.data import canonical_cds as ccds

# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions


# cds = cds.set_index(["seqname", "transcript_id", "exon_id"])
# cds["pos"] = cds.apply(lambda x: list(range(x["start"], x["end"] + 1)), axis=1)
# cds = cds.explode("pos")
# cds["cds_len"] = cds.groupby(level="transcript_id")["pos"].transform("count")
# cds["exon_len"] = cds.groupby(level="exon_id")["pos"].transform("count")


# # ## NMD annotations for + transcripts

# fwd = cds[cds["strand"] == "+"].copy()

# # ### Start proximal (distance from start codon)

# fwd["start_distance"] = (
#     fwd.groupby(level="transcript_id")["pos"].rank(ascending=True).astype(int)
# )

# # ### Long exons (distance upstream from splice junction)

# fwd["splice_donor_distance"] = (fwd["exon_end"] - fwd["pos"]) + 1

# # Where the CDS is in the last exon, there is no downstream splice donor site.
# # Here I will drop the "splice_donor_distance" annotation:
# fwd.loc[fwd["exon_number"] == fwd["exon_count"], "splice_donor_distance"] = np.nan

# # ### Last exon

# fwd["last_exon"] = np.where(fwd["exon_count"] == fwd["exon_number"], 1, 0)

# # ### 50nt rule

# fwd["fifty_nt_rule"] = np.where(
#     (fwd["exon_number"] == fwd["exon_count"] - 1)
#     & ((fwd["exon_end"] - fwd["pos"]) + 1 <= 50),
#     1,
#     0,
# )

# # ## NMD annotations for - transcripts

# rev = cds[cds["strand"] == "-"].copy()

# # ### Start proximal (distance from start codon)

# rev["start_distance"] = (
#     rev.groupby(level="transcript_id")["pos"].rank(ascending=False).astype(int)
# )

# # ### Long exons (distance upstream from splice junction)

# rev["splice_donor_distance"] = (rev["pos"] - rev["exon_start"]) + 1

# # Where the CDS is in the last exon (i.e. contiguous with the 3' UTR)
# # we should drop the "splice_donor_distance" annotation:
# rev.loc[rev["exon_number"] == rev["exon_count"], "splice_donor_distance"] = np.nan

# # ### Last exon

# rev["last_exon"] = np.where(rev["exon_count"] == rev["exon_number"], 1, 0)

# # ### 50nt rule

# rev["fifty_nt_rule"] = np.where(
#     (rev["exon_number"] == rev["exon_count"] - 1)
#     & ((rev["pos"] - rev["exon_start"]) + 1 <= 50),
#     1,
#     0,
# )

# # ## Merge fwd and rev annotations

# df = pd.concat([fwd, rev])
# df = df.reset_index()

# # ## Unify NMD annotations

# # Describe sites with overlapping NMD annotations
# a = pd.Series(np.where(df["start_distance"] <= 150, "start_proximal,", ""))
# b = pd.Series(np.where(df["splice_donor_distance"] > 400, "long_exon,", ""))
# c = pd.Series(np.where(df["fifty_nt_rule"] == 1, "fifty_nt,", ""))
# d = pd.Series(np.where(df["last_exon"] == 1, "last_exon,", ""))

# df["nmd"] = pd.Series(["".join([w, x, y, z]) for w, x, y, z in zip(a, b, c, d)])

# # Sites with no NMD-escape annotation are NMD targets
# df["nmd"] = df["nmd"].replace("", "nmd_target")

# # Create a definitive NMD annotation:
# ## "last_exon" and "fifty_nt" annotations are merged into "distal_nmd"
# ## Overlapping annotations are removed, with this priority:
# ## "start_proximal" > "distal_nmd" > "long_exon"
# df["nmd_definitive"] = df["nmd"].copy()
# df.loc[
#     df["nmd_definitive"].str.contains("start_proximal"), "nmd_definitive"
# ] = "start_proximal"
# df.loc[
#     (df["nmd_definitive"].str.contains("fifty_nt"))
#     | (df["nmd_definitive"].str.contains("last_exon")),
#     "nmd_definitive",
# ] = "distal_nmd"
# df["nmd_definitive"] = df["nmd_definitive"].replace("long_exon,", "long_exon")

# print("The number of positions per NMD annotation, genome-wide:")
# print(df.nmd_definitve.value_counts())


# df = (df[["seqname", "pos", "transcript_id", "nmd", "nmd_definitive"]]).rename(
#     columns={"seqname": "chr"}
# )
# df.to_csv("../outputs/nmd_annotations.tsv", sep="\t", index=False)
