#!/bin/bash

# Annotate ClinVar variants with VEP
# Run in the vep conda environment

vep \
    --input_file data/interim/clinvar_variants_selected.vcf \
    --output_file data/interim/clinvar_variants_vep.tsv \
    --species homo_sapiens \
    --assembly GRCh38 \
    --offline \
    --cache \
    --dir_cache /mnt/bmh01-rds/Ellingford_gene/.vep \
    --fork 8 \
    --buffer_size 50000 \
    --force_overwrite \
    --format vcf \
    --tab \
    --no_stats \
    --minimal \
    --coding_only \
    --show_ref_allele \
    --symbol \
    --fields "Location,REF_ALLELE,Allele,Consequence,Feature,Uploaded_variation,SYMBOL"
