#!/bin/bash

# Annotate SNVs with VEP

# Load conda environment for vep
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate
conda activate vep

# Run VEP
vep \
    --input_file data/interim/vep_all_snvs/in_$1.vcf \
    --output_file data/interim/vep_all_snvs/out_$1.tsv \
    --species homo_sapiens \
    --assembly GRCh38 \
    --offline \
    --cache \
    --dir_cache ~/.vep \
    --fork $(nproc) \
    --buffer_size 50000 \
    --force_overwrite \
    --format vcf \
    --tab \
    --no_stats \
    --minimal \
    --allele_number \
    --show_ref_allele \
    --coding_only \
    --fields "Location,REF_ALLELE,Allele,Consequence,Feature"

    # --per_gene \
    # --pick_order canonical \
