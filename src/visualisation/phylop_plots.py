"""Docstring."""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from src import constants as C
from src import visualisation as vis

# Styles and palettes
plt.style.use([C.STYLE_DEFAULT, C.COLOR_REGIONS])

# Functions
def horizontal_bars(
    values,
    ax=None,
    **kwargs,
):
    
    kwargs.setdefault("tick_label", values.index)
    kwargs.setdefault("color", sns.color_palette()[::-1])
    kwargs.setdefault("ecolor", [vis.adjust_lightness(c, 0.8) for c in kwargs["color"]])

    if not ax:
        ax = plt.gca()

    n = len(values)  # Number of bars
    y = np.arange(n)  # Y ticks

    bars = ax.barh(y=y / n, width=values, height=1 / (n + 1), **kwargs)

    return bars
