# Count Pfam domains overlapping NMD regions.
# Where an NMD region overlaps more than one domain, keep only the region /
# domain pair with the greatest overlap.
# Consider overlapping domains as those spanning >50% of the NMD region.

# TODO Get the total span of each NMD region
# TODO For each NMD region,

# %%
import pandas as pd

import src

FILE_IN = "data/interim/nmd_regions_pfam_intersect.tsv"
FILE_OUT = "data/final/pfam_domains_in_nmd_regions.tsv"

# %%
def keep_autosomal_loci(df):
    return df.loc[lambda x: ~x["chrom"].isin(["chrX", "chrY", "chrM"])]


def get_largest_domain_size(df):
    """Find the size of the largest overlapping domain for each NMD region."""

    return df.assign(
        region_domain_overlap=lambda x: x.groupby(
            ["enst", "region", "domain"], dropna=False
        )["overlap"]
        .transform("sum")
        .astype(int),
        greatest_domain_overlap=lambda x: x.groupby(["enst", "region"])[
            "region_domain_overlap"
        ].transform("max"),
    ).check.nrows(check_name="Rows after getting largest domain size")


def get_largest_domain_name(df):
    return (
        df.reset_index(drop=True)
        .assign(
            greatest_overlap_idx=lambda x: x.groupby(["enst", "region"])[
                "region_domain_overlap"
            ].transform(lambda g: g.idxmax()),
            domain_with_greatest_overlap=lambda x: x.loc[
                x["greatest_overlap_idx"], "domain"
            ].values,
        )
        .drop("greatest_overlap_idx", axis=1)
        .check.nrows(check_name="Rows after getting largest domain name")
    )


def get_region_sizes(df):
    return (
        df.drop_duplicates(["enst", "region", "start", "end"])
        .check.nrows(check_name="Rows after dropping duplicate exon/region pairs")
        .assign(
            exon_region_size=lambda x: x["end"] - x["start"],
            region_size=lambda x: x.groupby(["enst", "region"])[
                "exon_region_size"
            ].transform("sum"),
        )
    )


def choose_one_domain_per_region(df):
    """Keep regions without a constraint annotation."""
    return (
        df.groupby(["enst", "region", "constraint"], dropna=False)
        .agg(
            domain=("domain_with_greatest_overlap", "first"),
            overlap=("greatest_domain_overlap", "first"),
            region_size=("region_size", "first"),
        )
        .reset_index()
        .check.nrows(check_name="Rows after choosing one domain per region")
    )


def get_proportion_overlap(df):
    return df.assign(
        proportion_overlap=lambda x: x["overlap"] / x["region_size"]
    ).check.assert_data(lambda x: (x["proportion_overlap"] <= 1).all())


def main():
    return (
        pd.read_csv(
            FILE_IN,
            sep="\t",
            names=[
                "chrom",
                "start",
                "end",
                "enst",
                "ense",
                "region",
                "constraint",
                "domain",
                "overlap",
            ],
            na_values=".",
        )
        .pipe(keep_autosomal_loci)
        .pipe(get_largest_domain_size)
        .pipe(get_largest_domain_name)
        .pipe(get_region_sizes)
        .pipe(choose_one_domain_per_region)
        .pipe(get_proportion_overlap)
        .check.write(FILE_OUT, index=False)
        .check.head(5)
    )


if __name__ == "__main__":
    df = main()
