"""Visualisation tools for regional nonsense constraint paper."""

# Imports
from collections import namedtuple
from pathlib import Path

import colorsys
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import seaborn as sns

from src import constants as C
from src import setup_logger


# Logging
logger = setup_logger(Path(__file__).stem)


# Functions
def color_palette(style="default"):
    """Choose a color palette."""

    if not style in ["default", "regions", "maps"]:
        raise ValueError("style must be one of 'default', 'regions', or 'maps'.")

    if style == "default":
        labels = "blue green orange red light_blue pink grey black"
        plt.style.use(C.COLOR_VIBRANT)

    if style == "regions":
        labels = "transcript nmd_target start_proximal long_exon distal"
        plt.style.use(C.COLOR_REGIONS)

    if style == "maps":
        labels = C.MAPS_CONSEQUENCES
        plt.style.use(C.COLOR_MAPS)

    # Assign the color palette to a variable.
    # Colors can be selected by index or name (e.g. cp[0], cp.red)
    palette = namedtuple("color_palette", labels)
    palette = palette(*sns.color_palette().as_hex())

    return palette


def adjust_lightness(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))

    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def adjust_alpha(color, alpha):
    return mc.to_rgba(color, alpha)


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
