#!usr/bin/env bash

# This script annotates the ClinVar variants with VEP 105
# Run this script in the "vep" conda environment

vep \
    --input_file ../outputs/clinvar_variants_selected.vcf \
    --output_file ../outputs/clinvar_variants_vep.tsv \
    --species homo_sapiens \
    --assembly GRCh38 \
    --dir ~/.vep \
    --dir_cache ~/.vep \
    --offline \
    --cache \
    --fork 8 \
    --buffer_size 10000 \
    --force_overwrite \
    --format vcf \
    --tab \
    --no_stats \
    --pick \
    --pick_order mane,canonical \
    --coding_only