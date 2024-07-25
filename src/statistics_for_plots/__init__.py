"""Produce summary statistics for plots."""

import pandas as pd
from scipy import stats

_REGION_LABELS = {
    "transcript": "Full CDS",
    "full_cds": "Full CDS",
    "nmd_target": "NMD target",
    "start_proximal": "Start proximal",
    "long_exon": "Long exon",
    "distal_nmd": "Distal",
}


def sort_column(df, column="region", labels=_REGION_LABELS, **kwargs):
    """Rename and reorder the column of a dataframe."""

    kwargs.setdefault("ordered", True)

    assert all(
        i in labels for i in df[column]
    ), "Some values in the column are missing from the labels."

    df[column] = pd.Categorical(
        df[column], categories=labels, **kwargs
    ).rename_categories(labels)

    return df.sort_values(column)


def sort_index(df, labels=_REGION_LABELS, **kwargs):
    """Rename and reorder the index of a series or dataframe."""

    kwargs.setdefault("name", "region")

    assert all(
        i in labels.keys() for i in df.index
    ), "Some values in the index are missing from the labels."

    index = pd.Index(labels, **kwargs)

    return df.reindex(index).rename(index=labels)


def test_constrained_vs_unconstrained(df, **kwargs):
    """Welch's T-test of mean values in constrained vs unconstrained regions."""

    # Default arguments (can be overriden)
    kwargs.setdefault("equal_var", False)  # Welch's T test
    kwargs.setdefault("alternative", "two-sided")

    constrained = df.loc["constrained"]
    unconstrained = df.loc["unconstrained"]

    results = stats.ttest_ind_from_stats(
        constrained["mean"],
        constrained["std"],
        constrained["n"],
        unconstrained["mean"],
        unconstrained["std"],
        unconstrained["n"],
        **kwargs
    )

    return results
