#!/usr/bin/env bash
# Clean and index NMD region annotations

FILE_IN="data/interim/nmd_annotations.tsv"
TMP="data/final/nmd_annotations_simple.tsv"

tail -n +2 $FILE_IN \
| cut -f 1,2,3,20 \
| sort -k1,1V -k2,2n --stable -S 50% --parallel 8 \
> $TMP

# Index the file
bgzip --force $TMP
tabix --force --sequence 1 --begin 2 --end 2 "${TMP}.gz"