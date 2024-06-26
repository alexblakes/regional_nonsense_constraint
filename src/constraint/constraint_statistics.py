"""Identify regions and transcripts which are constrained for nonsense variants."""

import logging
from pathlib import Path

import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import fdrcorrection as fdr
from tqdm import tqdm

import src
from src import constants as C

_FILE_IN = "data/final/expected_variants_all_regions.tsv"
_FILE_OUT = "data/final/regional_constraint_stats.tsv"
_GNOMAD_V4_CONSTRAINT = "data/raw/gnomad.v4.0.constraint_metrics.tsv"
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_TRANSCRIPT = ["transcript"]
_REGIONS = ["distal_nmd", "nmd_target", "long_exon", "start_proximal"]
_CSQS = ["synonymous_variant", "missense_variant", "stop_gained"]

logger = logging.getLogger(__name__)


def per_row_binom_test(row):
    """Perform the binomial test per row."""

    try:
        result = binomtest(row["n_obs"], row["n_pos"], row["prop_exp"], alternative="less")
        p = result.pvalue

        # Get upper confidence interval of O/E value
        oe_ci_hi = (result.proportion_ci().high * row["n_pos"]) / row["n_exp"]

        return p, oe_ci_hi
    
    except: 
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

    # Read variant and mutability data
    df = pd.read_csv(
        _FILE_IN,
        sep="\t",
        # nrows=10000, #! Testing
    )

    # Get binomial p value and upper limit of O/E 95% CI
    logger.info("Running binomial tests per region and consequence type.")
    tqdm.pandas(desc="Per-row binomial tests.")
    binom_stat = df.progress_apply(per_row_binom_test, axis=1, result_type="expand").set_axis(
        ["p", "oe_ci_hi"], axis=1
    )

    df = pd.concat([df, binom_stat], axis=1)

    logger.info(f"Valid statistics: {df.p.count()}")
    logger.info(
        f"Available constraint statistics by region and csq:\n{df.groupby(['region', 'csq']).p.count()}"
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
    gnomad = get_gnomad_constraint(_GNOMAD_V4_CONSTRAINT)

    df = df.merge(gnomad, how="left")

    # Write to output
    df.to_csv(_FILE_OUT, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()
