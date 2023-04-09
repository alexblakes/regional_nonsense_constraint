#!/usr/bin/env bash

# This script downloads and combines the Hail coverage data at all VCF sites

# Make a local directory for the data
rm -rf ../data/split_coverage/
mkdir ../data/split_coverage

# Download the split_coverage files
dx find data --path outputs/gnomad_coverage/split_coverage |\
tr -s ' ' | cut -d ' ' -f 6 | head -n 100 |\
xargs -P 20 -I % dx download -o data/split_coverage/ %

# Concatenate the individual files containing the coverage data
# The number of files is too great for cat, so passing to xargs is required
cat ../data/split_coverage/*.tsv | sort -n -t t -k 1,2 > ../data/coverage.tsv

# Upload the combined data
dx rm -f outputs/coverage.tsv
dx upload --destination outputs/ ../data/coverage.tsv