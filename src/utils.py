"""Utility functions used throughout the analysis."""

import pandas as pd

from src.constants import REGION_LABELS

_REGION_LABELS = {
    "full_cds": "Full CDS",
    "nmd_target": "NMD target",
    "start_proximal": "Start proximal",
    "long_exon": "Long exon",
    "distal_nmd": "Distal",
}


def sort_index(series, labels=_REGION_LABELS, **kwargs):
    """Rename and order by the index for a series or dataframe."""
    kwargs.setdefault("name", "region")
    assert all(k in series.index for k in labels.keys()), "Index values are missing from labels.keys()"
    index = pd.Index(labels.keys(), **kwargs)
    return series.reindex(index=index).rename(labels)
