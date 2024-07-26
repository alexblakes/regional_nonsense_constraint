"""Visualisation tools for regional nonsense constraint paper."""

import colorsys

import matplotlib.colors as mc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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
    kwargs.setdefault("fontsize", 8)
    kwargs.setdefault("fontweight", "bold")
    kwargs.setdefault("transform", ax.transAxes)
    kwargs.setdefault("va", "bottom")
    kwargs.setdefault("ha", "right")

    ax.text(
        x,
        y,
        s,
        **kwargs,
    )

    return ax


def vertical_bars(series, ax=None, **kwargs):
    """Vertical bar chart for values in a series.

    xticklabels are taken from the series' index.
    """

    kwargs.setdefault("color", sns.color_palette())

    if not ax:
        ax = plt.gca()

    ax.bar(x=series.index, height=series, **kwargs)

    return ax


def vertical_grouped_bars(data, ax=None, bar_grouping="acmg", **kwargs):
    assert type(data) == pd.Series, "Data must be a pandas series"
    assert data.index.nlevels == 2, "Data must have a multiindex with two levels"

    if not ax:
        ax = plt.gca()

    bar_labels = data.index.get_level_values(bar_grouping).unique()

    for i, label in enumerate(bar_labels):
        data_subset = data.xs(label, level=bar_grouping)
        n_clusters = len(data_subset)
        x_position = np.arange(n_clusters)
        n_bars = len(bar_labels)
        bar_width = 1 / (n_bars + 1)
        offset = bar_width * i

        bars = ax.bar(
            x=x_position + offset,
            height=data_subset,
            width=bar_width,
            label=label,
            **kwargs,
        )

        # Place xticks centrally below the bar grouping...
        ax.set_xticks(x_position + bar_width, labels=data_subset.index)
        # ...But turn off xticklabels by default. 
        # To turn on elsewhere, use: 
        #   ax.tick_params(labelbottom=True)
        ax.tick_params(labelbottom=False)

    return ax
