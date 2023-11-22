"""Identify regions and transcripts which are constrained for nonsense variants."""


# Imports
from pathlib import Path

import pandas as pd
import numpy as np
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.proportion import proportions_chisquare
from statsmodels.stats.multitest import fdrcorrection as fdr

from src import setup_logger
from src import constants as C


# Module constants
_TRANSCRIPT = ["transcript"]
_REGIONS = ["distal_nmd", "nmd_target", "long_exon", "start_proximal"]
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
        alternative="two-sided",
        prop_var=row["prop_exp"],
    )

    return stat


def per_row_chisquare(row):
    """Calculate chi-square goodness of fit test for proportions row-wise."""

    if row["n_exp"] >= 5:  # Condition for X2 goodness of fit

        # Get statistics
        chi2, p, table = proportions_chisquare(
            count=row["n_obs"],
            nobs=row["n_pos"],
            value=row["prop_exp"],
        )

        z = np.sqrt(chi2) # Equivalent to two-sided z test

        # Negative z scores where O/E < 1
        if row["oe"] < 1:
            z = z * -1

        return chi2, z, p

    pass


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
            "constraint_flags": "gnomad_flags",
        }
    )

    return df


def main():
    """Run as script."""

    # Read expected variant data
    df = pd.read_csv(C.EXPECTED_VARIANTS_ALL_REGIONS, sep="\t")

    # Get chi squared scores, z scores, and unadjusted p values
    logger.info(
        "Running Chi squared goodness of fit tests per region and consequence type."
    )

    chi2 = df.apply(per_row_chisquare, axis=1, result_type="expand").set_axis(
        ["chi2", "z", "p"], axis=1
    )

    df = pd.concat([df, chi2], axis=1)

    logger.info(f"Valid chi squared statistics: {df.chi2.count()}")
    logger.info(
        f"Available constraint statistics by region and csq:\n{df.groupby(['region', 'csq']).chi2.count()}"
    )

    # FDR adjustment.
    # Calculated separately for whole-transcripts and constrained regions, and for each
    # distinct consequence.
    logger.info("Getting FDR statistics.")

    fdr_results = pd.concat(
        [
            fdr_adjustment(df, region=r, csq=c)
            for c in _CSQS
            for r in [_TRANSCRIPT, _REGIONS]
        ]
    )

    df = df.join(fdr_results)
    logger.info(f"Valid FDR results: {df.fdr_p.count()}")

    # Merge with gnomAD constraint data
    gnomad = get_gnomad_constraint(C.GNOMAD_V4_CONSTRAINT)

    df = df.merge(gnomad, how="left")

    # Write to output
    df.to_csv(C.REGIONAL_CONSTRAINT_STATS, sep="\t", index=False)

    # return df  #! Testing


if __name__ == "__main__":
    main()
