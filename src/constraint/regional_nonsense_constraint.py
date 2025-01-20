"""Find constrained and unconstrained regions."""



import pandas as pd

import src

_FILE_IN = "data/final/regional_constraint_stats.tsv"
_FILE_OUT = "data/final/regional_nonsense_constraint.tsv"

_OE_CI_HI = 0.6
_P = 0.05
_OBS = 1

logger = src.logger


def get_synonymous_p_vals(df):
    """Add synonymous p values to nonsense constraint data."""

    # Get synonymous and nonsense variants
    syn = df[df.csq == "synonymous_variant"].copy().rename(columns={"p": "syn_p"})
    stop = df[df.csq == "stop_gained"].copy()

    # Annotate nonsense constraint data with synonymous z scores
    stop = stop.merge(syn[["enst", "region", "syn_p"]], how="inner", validate="1:1")

    return stop


def assign_constraint_label(df):
    """Label constrained and unconstrained regions."""

    # Mark constrained regions
    m1 = df["oe_ci_hi"] <= _OE_CI_HI
    m2 = df["syn_p"] >= _P
    m3 = df["fdr_p"] < _P

    df.loc[m1 & m2 & m3, "constraint"] = "constrained"

    # Mark unconstrained regions
    m4 = df["p"] >= _P
    m5 = df["n_obs"] >= _OBS

    df.loc[m4 & m5, "constraint"] = "unconstrained"

    # Mark "indeterminate" regions
    m6 = df["fdr_p"].notna()
    m7 = df["constraint"].isna()

    df.loc[m6 & m7, "constraint"] = "indeterminate"

    # Logging
    logger.info(f"Constraint value counts:\n{df.constraint.value_counts(dropna=False)}")
    logger.info(
        f"Constrained region value counts:\n"
        f"{df.groupby('region').constraint.value_counts(dropna=False)}"
    )

    return df


def main():
    """Run as script."""

    # Read the data
    df = (
        pd.read_csv(_FILE_IN, sep="\t")
        .pipe(get_synonymous_p_vals)
        .pipe(assign_constraint_label)
    )

    # Write to output
    df.to_csv(_FILE_OUT, sep="\t", index=False)

    return df  #! Testing


if __name__ == "__main__":
    src.add_log_handlers()
    main()
