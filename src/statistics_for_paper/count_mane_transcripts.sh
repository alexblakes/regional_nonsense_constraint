#!/usr/bin/env bash
set -euo pipefail

# Count the number of MANE Select transcripts meeting our criteria.

_GTF="data/raw/gencode.v39.annotation.gtf"

grep -v "^#" $_GTF \
| awk '$3 == "transcript"' \
| grep 'transcript_type "protein_coding"' \
| grep 'gene_type "protein_coding"' \
| grep 'tag "Ensembl_canonical"' \
| tee >(wc -l | xargs echo "Canonical transcripts:" >&2) \
| grep 'tag "MANE_Select"' \
| tee >(wc -l | xargs echo "MANE Select transcripts:") \
> /dev/null