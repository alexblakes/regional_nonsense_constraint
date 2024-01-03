"""Colour tools."""

# Imports
from collections import namedtuple
from pathlib import Path

import colorsys
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import seaborn as sns

from src import constants as C
from src import setup_logger

logger = setup_logger(Path(__file__).stem)


# Useful functions
def colour_palette():
    """Define the default colour palette."""

    # Choose style sheet(s)
    plt.style.use(C.DEFAULT_MPL)
    sns.color_palette()

    # Assign the color palette to a variable.
    # Individual colors can be selected by index (e.g. C[0])
    cp = namedtuple("color_palette", "blue green orange red light_blue pink grey black")
    cp = cp(*sns.color_palette().as_hex())

    return cp


def adjust_lightness(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def panel_label(ax, s, x=-0.05, y=1.05, **kwargs):
    ax.text(
        x,
        y,
        s,
        transform=ax.transAxes,
        va="bottom",
        ha="right",
        fontsize=8,
        fontweight="bold",
        **kwargs,
    )


# Examples
## Qualitative

# # "Vibrant"
# # This is the default palette in the default.mplstyle style sheet
# # NB the red and magenta colours are hard to distinguish
# sns.color_palette(
#     [
#         "#0077BB",
#         "#009988",
#         "#EE7733",
#         "#CC3311",
#         "#33BBEE",
#         "#EE3377",
#         "#BBBBBB",
#         "#000000",
#     ]
# )


# # "Bright"
# # An alternative colour scheme where colours are easier to distinguish
# # However, it lacks a good red
# sns.color_palette(
#     ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"]
# )


# # The seaborn colorblind palette
# sns.color_palette("colorblind")


## Sequential

# sns.color_palette("Reds", as_cmap=True)

# # Using the "_r" suffix reverses the direction
# sns.color_palette("Blues_r", as_cmap=True)

# # Can also be a discrete scale, if the as_cmap argument is False (default)
# sns.color_palette("Blues_r")


## Divergent

# sns.color_palette("RdBu", as_cmap=True)

# # The number of values in a discrete scale can be set with the n_colors argument
# sns.color_palette("RdBu_r", n_colors=9)
