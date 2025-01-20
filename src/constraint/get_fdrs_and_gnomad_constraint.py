"""Identify regions and transcripts which are constrained for nonsense variants."""



import pandas as pd
from statsmodels.stats import multitest

import src

_FILE_IN = "data/interim/oe_stats_regions_cov_20.tsv"
_FILE_OUT = "data/final/regional_constraint_stats.tsv"
_GNOMAD_V4_CONSTRAINT = "data/raw/gnomad.v4.0.constraint_metrics.tsv"

_TRANSCRIPT = ["transcript"]
_REGIONS = ["distal_nmd", "nmd_target", "long_exon", "start_proximal"]
_CSQS = ["synonymous_variant", "missense_variant", "stop_gained"]

logger = src.logger


def fdr_adjustment(df, region, csq):
    """Get FDR-adjusted P-values for a given region and variant consequence."""

    # Mask regions and consequences
    m1 = df.region.isin(region)
    m2 = df.csq == csq

    # Filter the dataframe and drop cases lacking constraint statistics
    df = df.loc[m1 & m2, ["region", "p"]].dropna().copy()

    # FDR adjustment
    df["fdr_p"] = multitest.fdrcorrection(pvals=df["p"])[1]

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

    logger.info("Reading regional O/E statistics.")
    df = pd.read_csv(
        _FILE_IN,
        sep="\t",
        # nrows=10000,  #! Testing
    )

    logger.info(f"Observations: {len(df)}")
    logger.info(f"Valid P values: {df.p.count()}")
    logger.info(
        f"Available P values by region and csq:\n"
        f"{df.groupby(['region', 'csq']).p.count()}"
    )

    logger.info(
        "Getting FDR statistics. FDRs are calculated separately per consequence / "
        "region."
    )

    fdr_results = pd.concat(
        [
            fdr_adjustment(df, region=r, csq=c)
            for c in _CSQS
            for r in [_TRANSCRIPT, _REGIONS]
        ]
    )

    df = df.join(fdr_results, validate="1:1")
    logger.info(f"Valid FDR P values: {df.fdr_p.count()}")

    # Merge with gnomAD constraint data
    gnomad = get_gnomad_constraint(_GNOMAD_V4_CONSTRAINT)

    df = df.merge(gnomad, how="left")

    logger.info("Writing to output")
    df.to_csv(_FILE_OUT, sep="\t", index=False)

    return df


if __name__ == "__main__":
    src.add_log_handlers()
    main()
