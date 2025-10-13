#!/usr/bin/env bash
# Make BED file of NMD regions
set -eu

source src/utils.sh

FILE_NMD_POSITIONS="data/interim/nmd_annotations.tsv"
FILE_OUT="data/interim/nmd_regions_pre_bed.tsv"

< $FILE_NMD_POSITIONS tail -n+2 \
| tee >(_log "Input lines: " "$(wc -l)") \
| awk -v FS="\t" -v OFS="\t" '{print $1, $5, $3, $4, $20, $2}' \
| sort -k1,1V -k2,2 -k3,3V -k4,4V -k5,5 --stable -S 100M --parallel 16 \
| bedtools groupby -g 1-5 -c 6 -o min,max \
| awk -v FS="\t" -v OFS="\t" '{$6=$6 -1}1' \
| tee >(_log "Output lines (bed regions): " "$(wc -l)") \
> $FILE_OUT
