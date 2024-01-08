"""Utility functions to tidy data for plots."""

# Imports
import pandas as pd

from src.visualisation import plotting_constants as PC

# Module constants


# Functions
def categorical_regions_column(series):
    return pd.Categorical(
        series,
        categories=PC.REGIONS,
        ordered=True,
    ).rename_categories(PC.REGION_LABELS)


def categorical_regions_index(index, name="region"):
    return pd.CategoricalIndex(
        index,
        categories=PC.REGIONS,
        ordered=True,
        name=name,
    ).rename_categories(PC.REGION_LABELS)


def sort_region_column(df, column="region"):
    df[column] = categorical_regions_column(df[column])
    return df.sort_values("region")


def sort_region_index(df):
    df.index = categorical_regions_index(df.index)
    return df.sort_index()
