"""Docstring."""

# Imports
import numpy as np
import matplotlib.pyplot as plt


# Functions
def vertical_bars(df, ax, height, ylabel, bar_label=True, **kwargs):
    n = len(df)  # Number of bars
    x = np.arange(n)  # X ticks
    h = df[height]  # Height of bars (y-axis value)

    b = ax.bar(x=x / n, height=h, width=1 / (n + 1), **kwargs)

    ax.set_ylabel(ylabel)

    if bar_label:
        ax.bar_label(b, fmt="%.2f")

    # Add invisible x ticks
    ax.set_xticks(ticks=x / n, labels=[])
    ax.tick_params(axis="x", length=0)

    return None


def xticks(labels=[], ax=None, **kwargs):
    if ax == None:
        ax = plt.gca()

    if not labels:
        ax.set_xticks([])

    if labels:
        ax.set_xticks(
            ticks=ax.get_xticks(),
            labels=labels,
            rotation=45,
            rotation_mode="anchor",
            ha="right",
            **kwargs,
        )

    return None
