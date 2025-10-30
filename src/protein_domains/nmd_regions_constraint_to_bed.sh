#!/usr/bin/env bash
# Convert NMD regions / constraint annotation to BED format
set -eu

FILE_IN="data/interim/nmd_regions_pre_bed_constraint.tsv"
FILE_OUT="data/interim/nmd_regions_constraint.bed"

< $FILE_IN tail -n +2 \
| awk -v OFS="\t" '{print $1, $6, $7, ".", ".", $2, $3, $4, $5, $8}' \
> $FILE_OUT