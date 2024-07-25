"""Count VEP-annotated ClinVar variants."""

from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{'.'.join(Path(__file__).with_suffix('.log').parts[-2:])}"
_CLINVAR = "data/interim/clinvar_variants_vep_tidy.tsv"
_HEADER = "chr pos enst ref alt csq region acmg".split()

logger = src.setup_logger(_LOGFILE)

df = pd.read_csv(_CLINVAR, sep="\t", header=None, names=_HEADER)

logger.info(f"Valid ClinVar variants: {len(df)}")
logger.info(f"Consequence value counts:\n{df.csq.value_counts()}")