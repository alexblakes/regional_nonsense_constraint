"""Split ACMG and constraint annotations by NMD region."""

import argparse
from pathlib import Path

import pandas as pd

import src

logger = src.logger


def get_output_path(input_path, text_extension):
    no_suffix = Path(input_path).with_suffix("")
    old_stem = no_suffix.stem
    new_stem = old_stem + f"_{text_extension}"
    return no_suffix.with_name(new_stem).with_suffix(".tsv.gz")


def main():
    """Run as script."""

    parser = argparse.ArgumentParser()
    parser.add_argument("file_in", type=str)
    args = parser.parse_args()

    df = src.read_data(args.file_in).groupby("region")

    for region, data in df:
        file_out = get_output_path(args.file_in, region)
        src.write_out(data, file_out)

    return df


if __name__ == "__main__":
    src.add_log_handlers()
    df = main()
