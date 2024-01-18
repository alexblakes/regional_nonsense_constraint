"""Create a KDE plots for Z score distributions."""

# Imports
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
import numpy as np


# Functions
def z_dist(df, ax=None, x="z", xmin=-5, xmax=5, **kwargs):
    if not ax:
        ax = plt.gca()

    # Create the plot
    ax = sns.kdeplot(data=df, ax=ax, x=x, gridsize=1000, fill=False, **kwargs)

    # Adjust axes
    ax.set_xlim(xmin, xmax)

    return None


def axis_labels(title, ax=None, xlabel="Z score", ylabel=None):
    if not ax:
        ax = plt.gca()

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    return None


def extract_kde_line(ax=None):
    if not ax:
        ax = plt.gca()

    # Extract KDE line data
    kdeline = ax.lines[0]
    xs = kdeline.get_xdata()
    ys = kdeline.get_ydata()

    return xs, ys


def fill_kde(ax=None, **kwargs):
    if not ax:
        ax = plt.gca()

    xs, ys = extract_kde_line(ax)

    # Add fill beneath KDE line
    ax.fill_between(xs, 0, ys, **kwargs)

    return None


def fdr_line(df, ax=None, color="black", label="FDR < 0.05"):
    if not ax:
        ax = plt.gca()

    significant_results = df["fdr_p"] < 0.05
    max_p = df[significant_results]["p"].max()
    sig_thresh = norm.isf(1 - max_p)

    # Find y intercept for FDR threshold
    xs, ys = extract_kde_line(ax)
    y_val = np.interp(sig_thresh, xs, ys)

    # Add line
    ax.vlines(
        x=sig_thresh,
        ymin=0,
        ymax=y_val,
        color=color,
        linestyle="--",
        label=label,
    )

    if label:
        ax.text(
            x=sig_thresh,
            y=0.01,
            ha="right",
            va="bottom",
            s=f"{label}",
            rotation=90,
        )

    return None


def add_count(df, ax=None):
    if not ax:
        ax = plt.gca()

    N = len(df)
    ax.text(x=1, y=0.9, s=f"N = {N:,}", ha="right", va="top", transform=ax.transAxes)

    return None
