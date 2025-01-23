"""Visualisation tools for regional nonsense constraint paper."""

import colorsys

import matplotlib.colors as mc
import matplotlib.pyplot as plt
from matplotlib import transforms
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


def add_significance_asterisk(xs, ys, ps, ax=None, x_adj=0, y_adj=0, **kwargs):
    kwargs.setdefault("marker", (6, 2))
    kwargs.setdefault("color", "black")
    kwargs.setdefault("linewidth", 0.5)
    kwargs.setdefault("zorder", 3)
    kwargs.setdefault("clip_on", False)

    if not ax:
        ax = plt.gca()

    trans = ax.transData + transforms.ScaledTranslation(
        x_adj / 72, y_adj / 72, plt.gcf().dpi_scale_trans
    )

    for x, y, p in zip(xs, ys, ps):
        if p:
            ax.scatter(x, y, transform=trans, **kwargs)


def vertical_bars(series, ax=None, **kwargs):
    """Vertical bar chart for values in a series.

    xticklabels are taken from the series' index.
    """

    kwargs.setdefault("color", sns.color_palette())
    kwargs.setdefault("ecolor", [adjust_lightness(c, 0.8) for c in sns.color_palette()])

    if not ax:
        ax = plt.gca()

    ax.bar(x=series.index, height=series, **kwargs)

    return ax


def vertical_grouped_bars(data, ax=None, bar_grouping="acmg", yerrs=None, **kwargs):
    assert type(data) == pd.Series, "Data must be a pandas series"
    assert data.index.nlevels == 2, "Data must have a multiindex with two levels"
    if not yerrs is None:
        assert all(yerrs.index == data.index), "xerrs must have the same index as data"

    ax = ax or plt.gca()

    bar_labels = data.index.get_level_values(bar_grouping).unique()

    for i, label in enumerate(bar_labels):
        data_subset = data.xs(label, level=bar_grouping)
        n_clusters = len(data_subset)
        x_position = np.arange(n_clusters)
        n_bars = len(bar_labels)
        bar_width = 1 / (n_bars + 1)
        offset = bar_width * i

        if not yerrs is None:
            kwargs.update(yerr = yerrs.xs(label, level=bar_grouping))

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

def vertical_grouped_bars_test(data, bar_grouping, ax=None, yerr_column=None, **kwargs):
    if not isinstance(data, (pd.DataFrame, pd.Series)):
        raise TypeError("Data must be Pandas DataFrame or Series")
    if data.index.nlevels != 2:
        raise ValueError("Data must have a multiindex with exactly two levels")

    ax = ax or plt.gca()

