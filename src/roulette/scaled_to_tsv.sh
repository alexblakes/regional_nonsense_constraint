#!/usr/bin/env bash
set -euo pipefail

# Tidy and index the scaled Roulette scores, ready for annotation

FILE_IN="data/interim/roulette/scaled_cds_combined.tsv.gz"
FILE_OUT="data/interim/roulette/scaled_cds_combined_sorted.tsv.gz"
HEADER="data/manual/header_lines_roulette_scaled.txt"

# Write header for VCF
cat << EOF > $HEADER
##INFO=<ID=MR,Number=1,Type=Float,Description="Original Roulette mutation rate">
##INFO=<ID=MR_scaled,Number=1,Type=Float,Description="Scaled Roulette mutation rate in gnomAD v4.1">
EOF

zcat $FILE_IN \
| sort -u -k1,1V -k2,2n -k4,4 --parallel 8 --buffer-size 20% \
| bgzip -f > $FILE_OUT

tabix -s1 -b2 -e2 -f $FILE_OUT