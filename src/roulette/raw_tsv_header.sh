#!/usr/bin/env bash
set -euo pipefail

# Write a header for the final TSV of raw Roulette scores.

HEADER="data/manual/header_lines_roulette_raw_tsv.txt"

cat << EOF > $HEADER
##INFO=<ID=PN,Number=1,Type=String,Description="Pentanucleotide context">
##INFO=<ID=MR,Number=1,Type=Float,Description="Roulette mutation rate estimate">
##INFO=<ID=AR,Number=1,Type=Float,Description="Adjusted Roulette mutation rate estimate">
##INFO=<ID=MG,Number=1,Type=Float,Description="gnomAD mutation rate estimate (Karczewski et al. 2020)">
##INFO=<ID=MC,Number=1,Type=Float,Description="Carlson mutation rate estimate (Carlson et al. 2018)">
##INFO=<ID=MR_final,Number=1,Type=Float,Description="Roulette mutation rate, using adjusted mutation rate where available.">
EOF