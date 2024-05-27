#!/usr/bin/env bash

source activate vep # Before set command
set -euo pipefail

# Annotate all possible SNVs with VEP

ALL_SNVS="data/interim/cds_all_possible_snvs.vcf.gz"
TMP="data/tmp"
FILE_OUT="${TMP}/$1.vcf.gz"

bcftools view -r $1 $ALL_SNVS \
| vep \
    --input_file STDIN \
    --output_file STDOUT \
    --species homo_sapiens \
    --assembly GRCh38 \
    --offline \
    --cache \
    --dir_cache /mnt/bmh01-rds/Ellingford_gene/.vep \
    --fork 4 \
    --buffer_size 50000 \
    --format vcf \
    --vcf \
    --no_stats \
    --canonical \
    --fields "CANONICAL,Feature,Consequence" \
| filter_vep \
    --filter "CANONICAL is YES" \
    --filter "Consequence regex synonymous|missense|stop_gained" \
    --only_matched \
| bcftools +split-vep -c - -d \
| bcftools annotate \
    --set-id "%Feature" \
| bcftools annotate -x "^INFO/Consequence" -o $FILE_OUT