# Annotate the regions TSV with constraint data

# %%
import numpy as np
import pandas as pd
import src

FILE_REGIONS = "data/interim/nmd_regions_pre_bed.tsv"
FILE_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
FILE_OUT = "data/interim/nmd_regions_pre_bed_constraint.tsv"
FILE_OUT_BED = "data/interim/nmd_regions_color.bed"


# %%
def main():
    constraint = pd.read_csv(
        FILE_CONSTRAINT, sep="\t", usecols=["enst", "region", "constraint"]
    ).check.head()

    return (
        pd.read_csv(
            FILE_REGIONS,
            sep="\t",
            header=None,
            names=["chrom", "strand", "enst", "ense", "region", "start", "end"],
        )
        .check.nrows(check_name="Bed regions")
        .check.print("Left merge with constraint annotations.")
        .merge(constraint, how="left", validate="many_to_one")
        .check.nrows(check_name="Bed rows after merging constraint")
        .check.nrows(lambda x: x.drop_duplicates("enst"), check_name="Unique ENST IDs")
        .check.value_counts("region", check_name="Region value counts:")
        .check.function(
            lambda x: x.groupby("region")["constraint"].value_counts(dropna=False),
            check_name="Constraint value counts by region:",
        )
        .check.write(FILE_OUT, index=False)
        .assign(
            score="0",
            thick_start=lambda x: x["start"],
            thick_end=lambda x: x["end"],
            item_rgb=lambda x: np.select(
                condlist=[
                    x["region"] == "nmd_target",
                    x["region"] == "start_proximal",
                    x["region"] == "long_exon",
                    x["region"] == "distal_nmd",
                ],
                choicelist=[
                    "173,61,61",
                    "197,201,252",
                    "117,127,250",
                    "25,41,247",
                ],
                default="1,1,1",
            ),
        )
        .loc[
            :,
            [
                "chrom",
                "start",
                "end",
                "enst",
                "score",
                "strand",
                "thick_start",
                "thick_end",
                "item_rgb",
            ],
        ]
        .check.head()
        .check.write(FILE_OUT_BED, format="tsv", index=False, header=False)
    )


if __name__ == "__main__":
    main()
