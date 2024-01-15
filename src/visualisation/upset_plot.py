"""Create an upset plot of intersecting constrained regions."""

# Imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import upsetplot

from src import constants as C
from src import visualisation as vis


# Module constants
_GREY = vis.color_palette("default").grey

_UPSET_ARGS = dict(
    max_subset_size=2500,
    show_counts=True,
    element_size=20,
    intersection_plot_elements=3,
    totals_plot_elements=3,
    facecolor=_GREY,
)


# Visualise the plot
def create_upset_data(df, indicators=C.REGION_LABELS):
    return upsetplot.from_indicators(indicators, data=df)


def upset_plot(upset_data, fig=None, **kwargs):
    if not fig:
        fig = plt.gcf()

    upsetplot.plot(upset_data, fig=fig, **kwargs)

    return None


def customise_totals(ax=None):
    if not ax:
        ax = plt.gca()

    ax.grid(False)
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    ax.set_title("Total constrained", loc="right", ha="right")

    return None


def customise_intersection(ax=None):
    if not ax:
        ax = plt.gca()

    ax.set_ylabel("Intersecting\nconstrained", va="bottom", loc="bottom")
    ax.grid(False)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
