#!/usr/bin/env bash

# Annotate ClinVar variants with VEP
# Run in the vep conda environment

FILE_IN="data/interim/clinvar_variants_annotated.vcf"
FILE_OUT="data/interim/clinvar_variants_vep.vcf"

vep \
    --input_file $FILE_IN \
    --species homo_sapiens \
    --assembly GRCh38 \
    --offline \
    --cache \
    --dir_cache /mnt/bmh01-rds/Ellingford_gene/.vep \
    --fork 8 \
    --buffer_size 50000 \
    --force_overwrite \
    --format vcf \
    --vcf \
    --no_stats \
    --minimal \
    --symbol \
    --canonical \
    --fields "Consequence,Feature,SYMBOL,CANONICAL" \
    --output_file $FILE_OUT 

# Flags for TSV output:
# --fields "Location,REF_ALLELE,Allele,Consequence,Feature,Uploaded_variation,SYMBOL"
# --output_file data/interim/clinvar_variants_vep.tsv \