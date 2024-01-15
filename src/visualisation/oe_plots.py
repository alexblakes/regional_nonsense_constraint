"""Create a scatter plot with least squares line for transcript-level O/E statistics."""

# Imports
import matplotlib.pyplot as plt
import seaborn as sns


# Functions
def oe(df, ax=None, x="n_exp", y="n_obs", drop_ttn=True, alpha=0.2):
    if not ax:
        ax = plt.gca()

    if drop_ttn:
        df = df[df.enst != "ENST00000589042"]

    sns.regplot(
        ax=ax,
        data=df,
        x=x,
        y=y,
        ci=None,
        scatter_kws=dict(alpha=alpha),
    )

    ax.set_xlim(ax.get_ylim())


def add_x_eq_y_line(ax=None, **kwargs):
    if not ax:
        ax = plt.gca()

    ax.axline((0, 0), (1, 1), color="grey", linestyle="--")

    return None


def axis_labels(title, xlabel="Expected", ylabel="Observed", ax=None):
    if not ax:
        ax = plt.gca()

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return None
