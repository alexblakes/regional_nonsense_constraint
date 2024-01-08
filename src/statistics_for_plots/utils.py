"""Utility functions to tidy data for plots."""

# Imports
import pandas as pd

# Module constants
_REGION_ORDER = [
    "transcript",
    "nmd_target",
    "start_proximal",
    "long_exon",
    "distal_nmd",
]
_REGION_LABELS = [
    "Whole transcript",
    "NMD target",
    "Start proximal",
    "Long exon",
    "Distal",
]


# Functions
def categorical_regions_column(series):
    return pd.Categorical(
        series,
        categories=_REGION_ORDER,
        ordered=True,
    ).rename_categories(_REGION_LABELS)


def categorical_regions_index(index, name="region"):
    return pd.CategoricalIndex(
        index,
        categories=_REGION_ORDER,
        ordered=True,
        name=name,
    ).rename_categories(_REGION_LABELS)


def sort_region_column(df, column="region"):
    df[column] = categorical_regions_column(df[column])
    return df.sort_values("region")


def sort_region_index(df):
    df.index = categorical_regions_index(df.index)
    return df.sort_index()
