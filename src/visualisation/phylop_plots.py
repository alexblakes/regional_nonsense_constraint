"""Docstring."""

# Imports
import numpy as np
import matplotlib.pyplot as plt


# Functions
def horizontal_bars(df, ax=None, width, xlabel=False, yticklabels=[], **kwargs):
    if not ax:
        ax = plt.gca()
    
    n = len(df)  # Number of bars
    y = np.arange(n)  # Y ticks
    w = df[width]  # Width of bars (x-axis value)

    ax.bar(y=y / n, width=w, height=1 / (n + 1), **kwargs)

    ax.set_yticks(ticks=y * 1 / n, labels=yticklabels)
    ax.tick_params(axis="y", length=0)

    if xlabel:
        ax.set_xlabel(xlabel)

    return None
