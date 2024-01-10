"""Docstring."""

# Imports
import matplotlib.pyplot as plt


# Module constants
_TICK_NAMES = [
    "synonymous",
    "missense",
    "nonsense",
    "nmd_target",
    "start_proximal",
    "long_exon",
    "distal",
]
_TICK_LABELS = [
    "Synonymous",
    "Missense",
    "Nonsense\n(whole transcript)",
    "Nonsense\n(NMD target)",
    "Nonsense\n(start proximal)",
    "Nonsense\n(long exon)",
    "Nonsense\n(distal)",
]


# Functions
def plot_maps(df, y, x, xerr, color, ax=None):
    # Get current axis if not specified
    if ax == None:
        ax = plt.gca()

    # Point range plot
    ax.scatter(y=y, x=x, c=color)
    ax.errorbar(y=y, x=x, xerr=xerr, c=color, linestyle="None")

    return None


def yticks(ticks=_TICK_NAMES, labels=_TICK_LABELS, **kwargs):
    # Get current axis if not specified
    if ax == None:
        ax = plt.gca()

    # Set yticks
    ax.set_yticks(ticks=ticks, labels=labels, **kwargs)

    return None
