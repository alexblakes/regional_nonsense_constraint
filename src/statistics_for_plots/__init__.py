"""Package docstring."""

# Imports
import pandas as pd
from scipy import stats

from src import constants as C


# Functions
def categorical_regions_column(series, categories=C.REGIONS, labels=C.REGION_LABELS):
    """Create a categorical column."""
    return pd.Categorical(
        series,
        categories=categories,
        ordered=True,
    ).rename_categories(labels)


def categorical_regions_index(
    index, name="region", categories=C.REGIONS, labels=C.REGION_LABELS
):
    """Create a categorical index."""
    return pd.CategoricalIndex(
        index,
        categories=categories,
        ordered=True,
        name=name,
    ).rename_categories(labels)


def sort_region_column(df, column="region", **kwargs):
    """Sort a DataFrame by a categorical column."""

    df[column] = categorical_regions_column(df[column], **kwargs)
    return df.sort_values(column)


def sort_region_index(df, **kwargs):
    """Sort a DataFrame by a categorical column."""

    df.index = categorical_regions_index(df.index, **kwargs)
    return df.sort_index()


def test_constrained_vs_unconstrained(df, **kwargs):
    """Welch's T-test of mean values in constrained vs unconstrained regions."""

    # Default arguments (can be overriden)
    kwargs.setdefault("equal_var", False) # Welch's T test
    kwargs.setdefault("alternative","two-sided")

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
