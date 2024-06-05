#!/usr/bin/env bash
set -euo pipefail

# Reformat gnomAD v4.1 exome coverage data for annotation
# The v4.1 coverage data is apparently identical to v4.0

FILE_IN="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v4.1/exomes/coverage/gnomad.exomes.v4.0.coverage.summary.tsv.bgz"
FILE_OUT="data/interim/gnomad_v4.1_median_coverage.tsv.gz"
HEADER="data/manual/header_lines_coverage.txt"

# Write header
cat << EOF > $HEADER
##INFO=<ID=Coverage,Number=1,Type=Integer,Description="Median coverage in gnomAD v4.1 exomes.">
EOF

zcat $FILE_IN \
| tail -n +2 \
| sed 's/\:/\t/' \
| cut -f 1,2,4 \
| bgzip -f > $FILE_OUT

tabix -s1 -b2 -e2 -f $FILE_OUT