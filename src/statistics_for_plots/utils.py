"""Utility functions to tidy data for plots."""

# Imports
import pandas as pd

from src import constants as C

# Module constants


# Functions
def categorical_regions_column(series, categories=C.REGIONS, labels=C.REGION_LABELS):
    return pd.Categorical(
        series,
        categories=categories,
        ordered=True,
    ).rename_categories(labels)


def sort_region_column(df, column="region", **kwargs):
    df[column] = categorical_regions_column(df[column], **kwargs)
    return df.sort_values(column)


def categorical_regions_index(
    index, name="region", categories=C.REGIONS, labels=C.REGION_LABELS
):
    return pd.CategoricalIndex(
        index,
        categories=categories,
        ordered=True,
        name=name,
    ).rename_categories(labels)


def sort_region_index(df, **kwargs):
    df.index = categorical_regions_index(df.index, **kwargs)
    return df.sort_index()
