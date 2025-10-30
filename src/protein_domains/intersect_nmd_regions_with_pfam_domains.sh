#!/usr/bin/env bash
# Intersect NMD regions with Pfam domains.
set -eu

source src/utils.sh #_log

FILE_PFAM="data/interim/pfam_domains.bed"
FILE_NMD_REGIONS="data/interim/nmd_regions_constraint.bed"
FILE_OUT="data/interim/nmd_regions_pfam_intersect.tsv"

bedtools intersect -a $FILE_NMD_REGIONS -b $FILE_PFAM -wao -s \
| tee >(_log "Output lines: " "$(wc -l)") \
| awk -v OFS="\t" '{print $1,$2,$3,$7,$8,$9,$10,$14,$17}' \
> $FILE_OUT