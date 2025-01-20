"""Tidy constraint annotations for an upset plot."""



import pandas as pd

import src

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/final/regional_nonsense_constraint.tsv"
_CANONICAL_TRANSCRIPTS = "data/final/transcript_list_all.txt"
_FILE_OUT = "data/statistics/constraint_upset_data.tsv"
_REGION_LABELS = {
    "transcript": "Full CDS",
    "nmd_target": "NMD target",
    "start_proximal": "Start proximal",
    "long_exon": "Long exon",
    "distal_nmd": "Distal",
}


def read_constraint_data(path):
    return pd.read_csv(
        path,
        sep="\t",
        usecols=["region", "enst", "constraint"],
    )


def read_canonical_transcripts(path=_CANONICAL_TRANSCRIPTS):
    return pd.read_csv(path, sep="\t")


def filter_canonical_transcripts(df):
    canonical_transcripts = read_canonical_transcripts().squeeze()

    return df.loc[lambda x: x["enst"].isin(canonical_transcripts)]


def find_constrained_regions(df):
    df = (
        df.replace("unconstrained", False)
        .replace("indeterminate", False)
        .replace("constrained", True)
        .fillna(False)
    )
    logger.info(f"Constrained regions:\n{df.sum()}")

    return df


def write_out(df, path):
    df.to_csv(path, sep="\t")
    return df


def main():
    """Run as script."""

    return (
        read_constraint_data(_FILE_IN)
        .pipe(filter_canonical_transcripts)
        .pivot(index="enst", columns="region", values="constraint")
        .rename(columns=_REGION_LABELS)
        .loc[:, _REGION_LABELS.values()]
        .pipe(find_constrained_regions)
        .pipe(write_out, _FILE_OUT)
    )


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()

else:
    logger = src.logger
