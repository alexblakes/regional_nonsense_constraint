import matplotlib.pyplot as plt
import seaborn as sns


def vertical_bars(series, ax=None, **kwargs):
    """Vertical bar chart for values in a series.

    x tick labels are taken from the series' index.
    """

    kwargs.setdefault("color", sns.color_palette())

    if not ax:
        ax = plt.gca()

    ax.bar(x=series.index, height=series, **kwargs)

    return ax
