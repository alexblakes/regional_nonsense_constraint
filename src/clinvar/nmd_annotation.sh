#!/usr/bin/env bash
# Clean and index NMD region annotations

FILE_IN="data/interim/nmd_annotations.tsv"
TMP="data/final/nmd_annotations_simple.tsv"
HEADER="data/manual/header_lines_nmd_annotation.txt"

tail -n +2 $FILE_IN \
| cut -f 1,2,3,20 \
| sort -k1,1V -k2,2n --stable -S 50% --parallel 8 \
> $TMP

bgzip --force $TMP
tabix --force --sequence 1 --begin 2 --end 2 "${TMP}.gz"

# Write header line for NMD annotation
echo '##INFO=<ID=REGION,Number=1,Type=String,Description="NMD region">' > $HEADER
