"""Docstring."""

# Imports
import numpy as np
import matplotlib.pyplot as plt

from src import constants as C
from src import visualisation as vis

# Styles and palettes
plt.style.use(C.STYLE_DEFAULT)
_PALETTE = vis.color_palette("regions")[::-1]


# Functions
def horizontal_bars(
    values,
    ax=None,
    **kwargs,
):
    
    kwargs.setdefault("tick_label", values.index)
    kwargs.setdefault("color", _PALETTE)
    kwargs.setdefault("ecolor", [vis.adjust_lightness(c, 0.8) for c in _PALETTE])

    if not ax:
        ax = plt.gca()

    n = len(values)  # Number of bars
    y = np.arange(n)  # Y ticks

    ax.barh(y=y / n, width=values, height=1 / (n + 1), **kwargs)

    return None
