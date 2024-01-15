"""Create a hexbin plot for nonsense Z vs LOEUF scores."""

# Imports
import matplotlib.pyplot as plt


# Functions
def hexbin(df, ax=None, x="z", y="loeuf", zmin=-10, zmax=10, **kwargs):
    if not ax:
        ax = plt.gca()

    df = df[df["z"].between(zmin, zmax)]

    ax.hexbin(df[x], df[y], mincnt=1, gridsize=40, **kwargs)

    return None


def axis_labels(ax=None, x="Nonsense Z score", y="LOEUF"):
    if not ax:
        ax = plt.gca()
    ax.set_xlabel(x)
    ax.set_ylabel(y)

    return None


def add_rho(rho, ax=None, x=0.95, y=0.05):
    if not ax:
        ax = plt.gca()

    ax.text(
        x=x,
        y=y,
        s=rf"rho = {rho}",
        ha="right",
        va="bottom",
        transform=ax.transAxes,
    )

    return None
