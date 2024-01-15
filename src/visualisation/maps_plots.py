"""Create a horizontal point-range plot for MAPS scores."""

# Imports
import matplotlib.pyplot as plt


# Functions
def plot_maps(y, x, xerr, color, ax=None):
    # Get current axis if not specified
    if ax == None:
        ax = plt.gca()

    # Point range plot
    ax.scatter(y=y, x=x, c=color)
    ax.errorbar(y=y, x=x, xerr=xerr, c=color, linestyle="None")

    return None