# Find the proportion of regions (grouped by constraint and region type) which
# overlap a protein domain of any size.

# %%
import pandas as pd
import numpy as np
import pandas_checks
from statsmodels.stats import proportion

from src import constants as C
import pandas_utils as pdu

FILE_IN = "data/final/pfam_domains_in_nmd_regions.tsv"
FILE_OUT = "data/statistics/proportion_nmd_regions_overlapping_pfam_domains.tsv"


def per_row_proportion_ci(row):
    count = row["count"]
    nobs = row["total"]

    return proportion.proportion_confint(count, nobs, alpha=0.05, method="normal")


def main():
    return (
        pd.read_csv(FILE_IN, sep="\t")
        .loc[lambda x: x["constraint"].isin(["constrained", "unconstrained"])]
        .groupby(["constraint", "region"])
        .agg(count=("any_overlap", "sum"), total=("any_overlap", "size"))
        .reset_index()
        .assign(prop=lambda x: x["count"] / x["total"])
        .pipe(
            pdu.assign_with_apply,
            per_row_proportion_ci,
            new_cols=["ci_lo", "ci_hi"],
        )
        .assign(
            err_hi=lambda x: x["ci_hi"] - x["prop"],
            err_lo=lambda x: x["prop"] - x["ci_lo"],
        )
        .replace({"region": C.NMD_REGIONS_DICT})
        .assign(
            region=lambda x: pd.Categorical(x["region"], C.REGION_LABELS, ordered=True),
            constraint=lambda x: pd.Categorical(
                x["constraint"], ["unconstrained", "constrained"], ordered=True
            ).rename_categories(str.capitalize),
        )
        .sort_values(["constraint", "region"])
        .check.write(FILE_OUT, index=False)
    )


if __name__ == "__main__":
    df = main()
