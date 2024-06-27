#!/usr/bin/env bash
set -euo pipefail

# Filter SNVs with coverage >= 20 and AC > 0

FILE_IN="data/interim/cds_all_possible_snvs_annotated.tsv.gz"
FILE_OUT="data/interim/maps_snvs_filtered.tsv"
COVERAGE=20
AC=0

# Write header to output file
printf "csq\tregion\tac\tmu\tmu_adj\tmu_scaled\tmu_gnomad" > $FILE_OUT

# Filter SNVs
#   awk filters by coverage and AC, and drops irrelevant columns
#   output is appended to the file

zcat $FILE_IN \
| tail -n +2 \
| awk -v OFS="\t" -v CVR=${COVERAGE} -v AC=${AC} '
    $8 >= CVR && $9 > AC {print $6, $7, $9, $13, $14, $15, $16}
' \
>> $FILE_OUT