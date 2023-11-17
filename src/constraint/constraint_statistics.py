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


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def main():
    """Run as script."""

    df = pd.read_csv(C.EXPECTED_VARIANTS_ALL_REGIONS, sep="\t")

    return df #! Testing    

if __name__ == "__main__":
    main()




# # ## Get constraint Z scores


# def per_row_ztest(row, statistic="z"):
#     if statistic == "z":
#         i = 0
#     elif statistic == "p":
#         i = 1

#     stat = proportions_ztest(
#         count=row["n_obs"],
#         nobs=row["n_pos"],
#         value=row["prop_exp"],
#         alternative="smaller",
#         prop_var=row["prop_exp"],
#     )[i]

#     return stat


# %%capture
# df["z"] = df.apply(per_row_ztest, statistic="z", axis=1)
# df["p"] = df.apply(per_row_ztest, statistic="p", axis=1)

# # Print summary data
# _ = df.groupby(["region", "csq"])["z"].count()
# print(f"Constraint statistics by region and consequence:\n{_}")


# # ## Get FDR-adjusted P-values
# # For nonsense variants only. Calculate separately for whole-transcripts and constrained regions


# def fdr_adjustment(df, region, csq):
#     """Get FDR-adjusted P-values for a given region and variant consequence."""
#     # Mask regions and consequences
#     m1 = df.region.isin(region)
#     m2 = df.csq == csq

#     # Filter the dataframe and drop cases without a P-value
#     _df = df.loc[m1 & m2, ["region", "p"]].dropna().copy()

#     # FDR adjustment
#     _df["fdr_p"] = fdr(pvals=_df["p"])[1]

#     return _df["fdr_p"]


# # FDR adjustment is done separately for transcripts and NMD regions
# # and for each distinct consequence
# r1 = ["transcript"]
# r2 = ["distal_nmd", "nmd_target", "long_exon"]
# csq = ["synonymous", "missense", "nonsense"]

# fdr_results = pd.concat([fdr_adjustment(df, region=r, csq=c) for c in csq for r in [r1, r2]])


# # Join FDR-adjusted p-values to the original dataframe
# df = df.join(fdr_results)


# # ## Merge with gnomAD constraint data


# gnomad = pd.read_csv(
#     "../data/supplementary_dataset_11_full_constraint_metrics.tsv",
#     sep="\t",
#     usecols=["transcript", "pLI", "oe_lof_upper"],
# ).rename(columns={"transcript": "enst", "pLI": "pli", "oe_lof_upper": "loeuf"})

# df = df.merge(gnomad, how="left")


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


