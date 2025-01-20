"""Group ENSG IDs into "families" of paralogs."""



import pandas as pd

import src

FILE_IN = "data/interim/ensembl_paralogs_filtered.tsv"
FILE_OUT = "data/interim/paralog_families.tsv"

logger = src.logger


def read_paralogs(path=FILE_IN):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["ensg_target", "ensg_paralog"],
    )

    # Logging hereafter
    _only_target = df.loc[
        ~df["ensg_target"].isin(df["ensg_paralog"]), "ensg_target"
    ].drop_duplicates()

    _only_paralog = df.loc[
        ~df["ensg_paralog"].isin(df["ensg_target"]), "ensg_paralog"
    ].drop_duplicates()

    logger.info(f"Paralog data shape: {df.shape}")
    logger.info(f"Unique targets: {df.ensg_target.nunique()}")
    logger.info(f"Unique paralogs: {df.ensg_paralog.nunique()}")
    logger.info(f"ENSG IDs unique to targets: {len(_only_target)}")
    logger.info(f"ENSG IDs unique to paralogs: {len(_only_paralog)}")

    return df


def get_paralog_families(df):
    uniq_ensg_ids = set(df.to_numpy().flatten())

    paralog_families = set()

    for ensg in uniq_ensg_ids:
        # Get paralogs for each gene
        ensg_is_target = df["ensg_target"] == ensg
        ensg_is_paralog = df["ensg_paralog"] == ensg
        paralogs = df.loc[ensg_is_target | ensg_is_paralog].to_numpy().flatten()

        paralogs = frozenset(paralogs)  # Drops duplicates, is hashable

        paralog_families.add(paralogs)  # Duplicate families are removed

    paralog_families = pd.Series(list(paralog_families))  # Allows piping

    # Sanity checks
    gene_in_one_family = [
        sum([g in f for f in paralog_families]) == 1 for g in uniq_ensg_ids
    ]
    assert all(
        gene_in_one_family
    ), "Each gene must be a member of exactly one paralog family"

    logger.info(f"Number of paralog families: {len(paralog_families)}")

    return paralog_families


def assign_family_ids(series):
    series.index += 1  # Index will become the ID. We don't want an ID of 0.

    df = (
        series.rename("ensg")
        .to_frame()
        .explode("ensg")
        .reset_index(drop=False, names="paralog_family_id")
        .loc[:, ["ensg", "paralog_family_id"]]
    )

    logger.info(
        f"Paralog family sizes:\n"
        f"{df.groupby('paralog_family_id').count().describe()}"
    )

    return df


def write_out(df, path=FILE_OUT):
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""

    families = (
        read_paralogs()
        .pipe(get_paralog_families)
        .pipe(assign_family_ids)
        .pipe(write_out)
    )

    return families


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
