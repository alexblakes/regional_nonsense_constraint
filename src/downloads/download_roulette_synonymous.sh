#!/usr/bin/env bash
set -euo pipefail

# Download synonymous variant background for roulette

DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs"
URL="http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/all_hq_synonymous_variants.tsv.gz"
FILE="${DIR}/$(basename ${URL})"

wget -P $DIR $URL
gunzip $FILE