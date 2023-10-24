#! usr/bin/env bash

# This script annotates the ClinVar variants with VEP 105
# Run this script in the "vep" conda environment

# head -n $1 data/interim/cds_all_possible_snvs.vcf | \
vep \
    --input_file data/interim/cds_all_possible_snvs.vcf \
    --output_file data/interim/cds_all_possible_snvs_vep.vcf \
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
    --per_gene \
    --pick_order canonical \
    --coding_only \
    --fields "Location,REF_ALLELE,Allele,Consequence,Feature"

echo "All done."