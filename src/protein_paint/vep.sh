#!/usr/bin/env bash

# Annotate ClinVar variants with VEP
# Run in the vep conda environment

GENE=$1
FILE_IN=$2
FILE_OUT="data/interim/clinvar_variants_vep_protein_paint_${GENE}.vcf"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

vep \
    --input_file $FILE_IN \
    --species homo_sapiens \
    --assembly GRCh38 \
    --fasta $FASTA \
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
    --hgvs \
    --fields "Consequence,Feature,SYMBOL,CANONICAL,HGVSp" \
    --output_file $FILE_OUT 

# Flags for TSV output:
# --fields "Location,REF_ALLELE,Allele,Consequence,Feature,Uploaded_variation,SYMBOL"
# --output_file data/interim/clinvar_variants_vep.tsv \