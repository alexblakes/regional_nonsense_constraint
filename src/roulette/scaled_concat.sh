#!/usr/bin/env bash
set -euo pipefail

# Concatenate the scaled Roulette mutation rate scores

export TMP="data/tmp"
FILE_OUT="data/interim/roulette/scaled_cds_combined.tsv.gz"
FILE_PATHS_FILTERED="${TMP}/roulette_filtered_paths.txt"
FILE_PATHS_SCALED="${TMP}/roulette_scaled_paths.txt"
export REGIONS="data/interim/gencode_v39_canonical_cds.bed"
ROULETTE_DIR="data/interim/roulette"

mkdir -p $TMP

# Write the paths of the raw Roulette scores to a text file.
find $ROULETTE_DIR -type f -name *.tsv.bgz | sort -k1,1V > $FILE_PATHS_SCALED

limit_to_regions() {
    # Only keep sites in CDS of canonical transcripts.

    local F_OUT="${TMP}/$(basename ${1})"
    tabix -R $REGIONS $1 | bgzip > $F_OUT
}

# Filter in parallel
export -f limit_to_regions
parallel --arg-file $FILE_PATHS_SCALED -j80% limit_to_regions {}

# Concatenate the files 
zcat $(find $TMP -type f -name *.tsv.bgz | sort -V) | bgzip -f > $FILE_OUT
# The positions are unsorted, so indexing is not possible yet.

# Clean up
rm -r $TMP