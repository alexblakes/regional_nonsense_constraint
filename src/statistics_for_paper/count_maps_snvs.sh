#!/usr/bin/env bash
set -euo pipefail

# Count SNVs used for the MAPS calculation

SNVS="data/interim/cds_all_possible_snvs_annotated.tsv.gz"

zcat $SNVS \
| tail -n +2 \
| awk '$9 > 0 && $8 >= 20' \
| cut -f 1,2,3,4 \
| uniq \
| wc -l | xargs echo "Unique SNVs qualifying for MAPS:"
