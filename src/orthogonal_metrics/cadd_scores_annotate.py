"""Annotate sites with CADD scores."""

import itertools


import numpy as np
import pandas as pd
from sklego.pandas_utils import log_step

import src
from src.constraint.observed_variants_counts_and_oe_stats import reorder_data
from src.orthogonal_metrics import merge_orthogonal_annotations as moa

_CADD = "data/interim/cadd_scores_coding.tsv"
_NMD = "data/interim/nmd_annotations.tsv"
_VEP = "data/interim/cds_all_possible_snvs_vep_tidy.tsv"
_CONSTRAINT = "data/final/regional_nonsense_constraint.tsv"
_FILE_OUT = "data/interim/cadd_scores_coding_annotated.tsv"

_DTYPES = {
    "chr": "category",
    "pos": np.int32,
    "ref": "category",
    "alt": "category",
    "cadd_phred": np.float16,
    "csq": "category",
    "enst": "category",
    "region": "category",
    "constraint": "category",
}

logger = src.logger


def read_cadd(path, **kwargs):
    """Read CADD data."""

    logger.info("Reading CADD data.")

    cadd = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names="chr pos ref alt cadd_raw cadd_phred".split(),
        usecols="chr pos ref alt cadd_phred".split(),
        # low_memory=False,
        dtype=_DTYPES,
        # nrows=10000,
        **kwargs,
    )

    return cadd


@log_step(print_fn=lambda x: logger.info(x), shape_delta=True)
def tidy_cadd(df):
    return (
        df.drop_duplicates()
        .assign(chr=lambda x: ["".join(["chr", str(c)]) for c in x["chr"]])
        .assign(chr=lambda x: pd.Categorical(x["chr"]))
    )


def read_vep(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "pos", "ref", "alt", "csq", "enst"],
        dtype={**_DTYPES},
        # nrows=10000,
    )


def set_dtypes(df):
    return df.astype({x: _DTYPES[x] for x in df.columns})


def harmonise_categories(dfs):
    get_cat_columns = lambda x: x.select_dtypes(include="category").columns
    cat_columns = [set(get_cat_columns(df)) for df in dfs]
    shared_cat_columns = set.intersection(*cat_columns)

    for col in shared_cat_columns:
        uc = pd.api.types.union_categoricals([df[col] for df in dfs])
        for df in dfs:
            df[col] = df[col].cat.set_categories(uc.categories)

    return dfs


@log_step(print_fn=lambda x: logger.info(x), shape_delta=True, display_args=False)
def merge_l_and_r(l, r):
    l, r = harmonise_categories([l, r])
    return l.merge(r, how="inner", validate="many_to_one")


def reorder_data(df):
    return df[["enst", "region", "csq", "constraint", "cadd_phred"]]


def write_out(df, path):
    df.to_csv(path, sep="\t", index=False)
    return df


def main():
    """Run as script."""

    # Read the data
    cadd = read_cadd(_CADD).pipe(tidy_cadd)
    vep = read_vep(_VEP)
    nmd = moa.read_nmd_annotations(_NMD).pipe(set_dtypes)
    constraint = moa.read_regional_nonsense_constraint(_CONSTRAINT).pipe(set_dtypes)

    df = (
        merge_l_and_r(vep, cadd)
        .pipe(merge_l_and_r, nmd)
        .pipe(merge_l_and_r, constraint)
        .pipe(reorder_data)
        .pipe(write_out, _FILE_OUT)
    )
    # # Write to output
    # logger.info("Writing to output.")
    # df.to_csv(_FILE_OUT, sep="\t", index=False)

    # logger.info("Done.")

    return df


if __name__ == "__main__":
    src.add_log_handlers()
    main()
