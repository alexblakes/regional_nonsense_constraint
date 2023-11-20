"""Constraint statistics."""


# Imports
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.multitest import fdrcorrection as fdr
from scipy import stats as _stats
from scipy.stats import spearmanr

from src import setup_logger
from src import constants as C


# Module constants
_TRANSCRIPT = ["transcript"]
_REGIONS = ["distal_nmd", "nmd_target", "long_exon"]
_CSQS = ["synonymous_variant", "missense_variant", "stop_gained"]


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def per_row_ztest(row):
    """Calculate one-sided z test row-wise"""

    stat = proportions_ztest(
        count=row["n_obs"],
        nobs=row["n_pos"],
        value=row["prop_exp"],
        alternative="smaller",
        prop_var=row["prop_exp"],
    )

    return stat


def fdr_adjustment(df, region, csq):
    """Get FDR-adjusted P-values for a given region and variant consequence."""

    # Mask regions and consequences
    m1 = df.region.isin(region)
    m2 = df.csq == csq

    # Filter the dataframe and drop cases lacking constraint statistics
    df = df.loc[m1 & m2, ["region", "p"]].dropna().copy()

    # FDR adjustment
    df["fdr_p"] = fdr(pvals=df["p"])[1]

    return df["fdr_p"]


def get_gnomad_constraint(path):
    """Read gnomAD v4 constraint data."""

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["transcript", "lof.pLI", "lof.oe_ci.upper", "constraint_flags"],
    ).rename(
        columns={
            "transcript": "enst",
            "lof.pLI": "pli",
            "lof.oe_ci.upper": "loeuf",
            "constraint_flags": "gnomad_constraint_flags",
        }
    )

    return df


def main():
    """Run as script."""

    # Read expected variant data
    df = pd.read_csv(C.EXPECTED_VARIANTS_ALL_REGIONS, sep="\t", nrows=5)

    # Get z scores and unadjusted p values
    zp = df.apply(per_row_ztest, axis=1, result_type="expand").set_axis(
        ["z", "p"], axis=1
    )

    df = pd.concat([df, zp], axis=1)

    logger.info(
        f"Available constraint statistics by region and csq:\n{df.groupby(['region', 'csq']).z.count()}"
    )

    # FDR adjustment.
    # Calculated separately for whole-transcripts and constrained regions, and for each
    # distinct consequence.
    fdr_results = pd.concat(
        [
            fdr_adjustment(df, region=r, csq=c)
            for c in _CSQS
            for r in [_TRANSCRIPT, _REGIONS]
        ]
    )

    df = df.join(fdr_results)

    # Merge with gnomAD constraint data
    gnomad = get_gnomad_constraint(C.GNOMAD_V4_CONSTRAINT)
    
    df = df.merge(gnomad, how="left")

    return df  #! Testing


if __name__ == "__main__":
    main()


# # ## Save to output


# df.to_csv("../outputs/expected_variants_all_regions_stats.tsv", sep="\t", index=False)


# # ---
# # ## Statistics


# # Re-load the data
# df = pd.read_csv("../outputs/expected_variants_all_regions_stats.tsv", sep="\t")


# # ### Exclude regions where the synonymous Z-score is < -1


# m1 = df.csq == "synonymous"
# m2 = df.z >= -1
# m = df[m1 & m2][["region", "enst"]]
# df2 = df.merge(m, how="inner")


# # ### Count regions in p-value bins


# def get_p_stats_region(df, region, csq="nonsense", bins=[0, 0.001, 0.01, 0.05, 1]):
#     m1 = df.region == region
#     m2 = df.csq == csq

#     df = df[m1 & m2][["n_obs", "p", "fdr_p"]]

#     for p in ["p", "fdr_p"]:
#         b = f"{p}_bin"
#         df[b] = pd.cut(df[p], bins=bins)
#         g = df.groupby(b)
#         stats = g.agg({f"{p}": "count"})
#         stats["none_observed"] = g["n_obs"].apply(lambda x: (x == 0).sum())
#         print(f"{csq}, {region}")
#         print(f"{stats}\n")


# for r in df.region.unique():
#     get_p_stats_region(df2, region=r)


# # ### Spearman rank Z vs LOEUF


# m1 = df["region"] == "transcript"
# m2 = df["csq"] == "nonsense"

# z = df[m1 & m2]["z"]
# loeuf = df[m1 & m2]["loeuf"]

# print(spearmanr(z, loeuf, nan_policy="omit", alternative="two-sided"))
# print("\nSee below for smallest possible P-value ('tiny')")
# print(np.finfo(np.float64))
