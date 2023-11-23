"""Docstring."""

# Imports
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as _stats
import statsmodels.formula.api as smf
from statsmodels.stats.proportion import proportions_ztest

from src import setup_logger
from src import constants as C

# Module constants
_VCF_HEADER = ["chr", "pos", "id", "ref", "alt", "qual", "filter", "info"]
_DATATYPES = defaultdict(lambda: "str")
_DATATYPES.update({"pos": np.int32, "ac": np.int32, "an": np.int32})
_COVERAGE = 20


# Logging
logger = setup_logger(Path(__file__).stem)

# Functions
def get_variant_annotations(path):
    """Get merged annotations for all variants."""

    df = (
        pd.read_csv(
            path,
            sep="\t",
            # nrows=1000000,
        )
        .rename(columns={"nmd": "region"})
        .query("obs == True")
        .query(f"median_coverage >= {_COVERAGE}")
    )

    # Sanity checks
    assert all(x > 0 for x in df.ac), "0s or NaNs present in df.ac"

    # Logging
    logger.info(f"Observed variants at {_COVERAGE}x coverage: {len(df)}")
    logger.info(f"Missing values:\n{df.isna().sum()}")
    
    return df


def get_valid_synonymous_variants(df):
    """Get valid synonymous variants for MAPS model."""

    logger.info("Filtering for valid synonymous variants.")

    syn = df[df["csq"] == "synonymous_variant"]

    logger.info(f"Synonymous variants before filtering: {len(syn)}")

    # Mask contexts in which a synonymous variant is generally not possible.
    # Synonymous variants in these contexts can only occur at exon-intron junctions.
    m1 = (syn.tri == "AGT") & ((syn.alt == "C") | (syn.alt == "T"))
    m2 = (syn.tri == "AAT") & ((syn.alt == "C") | (syn.alt == "T"))
    m3 = (syn.tri == "ACT") & ((syn.alt == "G") | (syn.alt == "A"))
    m4 = (syn.tri == "ATT") & ((syn.alt == "G") | (syn.alt == "A"))

    syn = syn[~(m1 | m2 | m3 | m4)]

    logger.info(f"Synonymous variants after filtering: {len(syn)}")

    return syn


def get_ps(dfg):
    """Get mean mutability and proportion of singletons."""

    mu = dfg["mu"].mean()
    ns = dfg["ac"].apply(lambda x: (x == 1).sum()).rename("n_singletons")
    no = dfg["ac"].count().rename("n_obs")

    ps = pd.concat([mu, ns, no], axis=1).reset_index()

    ps["ps"] = ps["n_singletons"] / ps["n_obs"]

    return ps


def main():
    """Run as script."""

    df = (
        get_variant_annotations(C.ALL_VARIANTS_MERGED_ANNOTATIONS)
        .pipe(get_valid_synonymous_variants)
        .groupby(["tri", "ref", "alt", "variant_type", "lvl"])
        .pipe(get_ps)
    )

    # Write to output
    df.to_csv(C.PS_SYN_CONTEXT, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    main()


# ps.to_csv("../outputs/proportion_singletons_synonymous_by_context.tsv", sep="\t", index=False)


# # ### Proportion of singletons by csq


# m1 = df["ac"] != 0
# m2 = df["variant_type"] != "CpG"
# m3 = df["csq"] == "nonsense"


# # #### Synonymous, missense, and nonsense


# # CpG only
# ps_csq_cpg = get_ps(df[m1 & ~m2].groupby("csq"))

# # Exclude CpG
# ps_csq_no_cpg = get_ps(df[m1 & m2].groupby("csq"))


# # #### Nonsense, by NMD region


# # CpG only
# ps_region_cpg = get_ps(df[m1 & m3 & ~m2].groupby("nmd"))

# # Exclude CpG
# ps_region_no_cpg = get_ps(df[m1 & m2 & m3].groupby("nmd"))


# # ### Combine CSQ and Region results


# # Reformat column names
# ps_region_cpg = ps_region_cpg.rename(columns={"nmd":"csq"})
# ps_region_no_cpg = ps_region_no_cpg.rename(columns={"nmd":"csq"})

# for a, b, c in zip([ps_csq_cpg, ps_csq_no_cpg], [ps_region_cpg, ps_region_no_cpg], ["cpg","no_cpg"]):
#     ps = pd.concat([a,b])
#     ps.to_csv(f"../outputs/proportion_singletons_by_csq_and_region_{c}.tsv", sep="\t", index=False)
#     print(f"{c}\n{ps}\n\n")
