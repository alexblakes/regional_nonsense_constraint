#!/usr/bin/env bash
set -euo pipefail

# Drop the header line from the scaled Roulette score TSVs

export SCALED_DIR="data/interim/roulette"
TMP="data/tmp_scaled"
FILE_PATHS="${TMP}/scaled_rate_paths.txt"

# Create temporary directory
mkdir -p $TMP

# Write file paths to a text file
find $SCALED_DIR -type f -name *_all.tsv.gz | sort -V > $FILE_PATHS

remove_header() {
    # Remove header
    # Add "chr" prefix
    
    FILE_IN=$1
    FILE_OUT="${SCALED_DIR}/$(basename -s .gz ${1}).bgz"

    zcat $FILE_IN \
    | tail -n +2 \
    | sed 's/^/chr/' \
    | bgzip > $FILE_OUT    
    tabix -s1 -e2 -b2 -f $FILE_OUT
}

# Process files in parallel
export -f remove_header
parallel --arg-file $FILE_PATHS -j80% remove_header {}

# Clean up
rm $TMP