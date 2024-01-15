"""Docstring."""

# Imports
import matplotlib.pyplot as plt

from src import constants as C

# Module constants



# Functions
def plot_maps(y, x, xerr, color, ax=None):
    # Get current axis if not specified
    if ax == None:
        ax = plt.gca()

    # Point range plot
    ax.scatter(y=y, x=x, c=color)
    ax.errorbar(y=y, x=x, xerr=xerr, c=color, linestyle="None")

    return None


def yticks(ticks=C.MAPS_LABELS, labels=C.MAPS_LABELS, **kwargs):
    # Get current axis if not specified
    if ax == None:
        ax = plt.gca()

    # Set yticks
    ax.set_yticks(ticks=ticks, labels=labels, **kwargs)

    return None
