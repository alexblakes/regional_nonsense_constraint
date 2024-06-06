#!/usr/bin/env bash

source activate vep # Before set command
set -euo pipefail

# Annotate all possible SNVs with VEP
# Sending the VCF to a compressed output leads to a bgzf_read error with bcftools concat.
# Therefore, the output is not compressed.

ALL_SNVS="data/interim/cds_all_possible_snvs.vcf.gz"
TMP="data/tmp"
FILE_OUT="${TMP}/$1.vcf"

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
    --buffer_size 5000 \
    --format vcf \
    --vcf \
    --no_stats \
    --canonical \
    --fields "CANONICAL,Feature,Consequence" \
    --warning_file "data/logs/vep_warn_${1}.txt" \
    --allow_non_variant \
    --dont_skip \
| filter_vep \
    --filter "CANONICAL is YES" \
    --filter "Consequence regex synonymous|missense|stop_gained" \
    --only_matched \
| bcftools +split-vep -c - -d -Ou \
| bcftools annotate --set-id "%Feature" -Ou \
| bcftools annotate -x "^INFO/Consequence" -o $FILE_OUT