#!/usr/bin/env bash
set -euo pipefail

# Index the bgzipped Roulette files

ROULETTE_DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs"
TMP="data/tmp"
FILE_PATHS="${TMP}/roulette_raw_paths.txt"

# Make a temporary directory
rm -r $TMP
mkdir -p $TMP

# Write file paths to a text file
find $ROULETTE_DIR -type f -name *.vcf.bgz | sort -k1,1V > $FILE_PATHS

# Index the files
parallel --arg-file $FILE_PATHS tabix -f {}

# Clean up
rm -r $TMP