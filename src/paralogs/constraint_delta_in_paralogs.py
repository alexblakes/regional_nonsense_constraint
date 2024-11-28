"""Find paralogous genes with strongly different constraint."""

import logging

import pandas as pd

import src

FILE_PARALOG_FAMILIES = "data/interim/paralog_families.tsv"
FILE_CONSTRAINT = "data/final/regional_nonsense_constraint_for_excel.tsv"
FILE_OUT = "data/interim/paralogs_delta_constraint.tsv"

logger = logging.getLogger(__name__)


def read_paralogs(path=FILE_PARALOG_FAMILIES):
    df = pd.read_csv(path, sep="\t")

    logger.info(
        f"Data shape: {df.shape}\n"
        f"Unique ENSGs: {df.ensg.nunique()}\n"
        f"Unique paralog families: {df.paralog_family_id.nunique()}"
    )

    return df


def read_constraint(path=FILE_CONSTRAINT):
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["symbol", "ensg", "enst", "region", "oe_ci_hi", "constraint"],
    )

    logger.info(
        f"Data shape: {df.shape}\n"
        f"Unique symbols {df.symbol.nunique()}\n"
        f"Unique ENSGs: {df.ensg.nunique()}\n"
        f"Unique ENSTs: {df.enst.nunique()}"
    )

    return df


def merge_constraint_paralogs(left, right):
    df = left.merge(right, how="inner", validate="many_to_one")

    logger.info(
        f"Data shape: {df.shape}\n"
        f"Unique ENSGs: {df.ensg.nunique()}\n"
        f"Unique paralog_families: {df.paralog_family_id.nunique()}"
    )
    return df


def filter_constrained(df):
    """Filter for paralog families in which at least one member is constrained in a
    given region.
    """

    df = df.groupby(["paralog_family_id", "region"]).filter(
        lambda x: any(x["constraint"] == "constrained")
    )

    logger.info(
        f"Data shape: {df.shape}\n"
        f"Unique ENSGs: {df.ensg.nunique()}\n"
        f"Unique paralog_families: {df.paralog_family_id.nunique()}\n"
        f"Number of groups by paralog / region: {len(df.drop_duplicates(['paralog_family_id','region']))}"
    )

    return df


def get_delta_constraint(df):
    grouped = df.groupby(["paralog_family_id", "region"])

    df = grouped.agg(
        paralogs=("symbol", ",".join),
        oe95_min=("oe_ci_hi", "min"),
        oe95_max=("oe_ci_hi", "max"),
    ).assign(oe95_delta=lambda x: x.oe95_max - x.oe95_min)

    # Can't fit these into named aggregation
    most_constrained = grouped.apply(
        lambda x: x.loc[x["oe_ci_hi"].idxmin(), "symbol"]
    ).rename("most_constrained")
    least_constrained = grouped.apply(
        lambda x: x.loc[x["oe_ci_hi"].idxmax(), "symbol"]
    ).rename("least_constrained")

    _shape_df_pre_merge = df.shape  # For logging

    df = (
        pd.concat([df, most_constrained, least_constrained], join="inner", axis=1)
        .loc[
            :,
            [
                "paralogs",
                "most_constrained",
                "least_constrained",
                "oe95_min",
                "oe95_max",
                "oe95_delta",
            ],
        ]
        .reset_index()
    )

    assert _shape_df_pre_merge[0] == df.shape[0], "Should be a one-to-one merge."

    logger.info(
        f"Shape before merge: {_shape_df_pre_merge}\n" f"Shape after merge: {df.shape}"
    )

    return df


def write_out(df, path=FILE_OUT):
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""

    paralogs = read_paralogs()
    constraint = read_constraint()

    df = (
        merge_constraint_paralogs(constraint, paralogs)
        .pipe(filter_constrained)
        .pipe(get_delta_constraint)
        .pipe(write_out)
    )

    return df


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
