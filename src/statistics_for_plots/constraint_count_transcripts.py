"""Count the number of transcripts for which regional constraints could be quantified."""



import pandas as pd

import src
from src import statistics_for_plots as sp

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_FILE_IN = "data/final/regional_nonsense_constraint.tsv"
_FILE_OUT = "data/statistics/constraint_count_transcripts.tsv"

logger = src.logger


def read_data(path):
    return pd.read_csv(
        path,
        sep="\t",
        usecols=["enst","region","n_obs","n_exp","oe","p","fdr_p","pli","loeuf"]
    )


def write_out(series, path):
    series.to_csv(path, sep="\t")
    return series


def main():
    """Run as script."""

    df = (
        read_data(_FILE_IN)
        # .pipe(write_out, _FILE_OUT)
    )

    return df


if __name__ == "__main__":
    src.add_log_handlers()
    main()
