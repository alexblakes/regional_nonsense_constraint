#!/bin/bash

# Concatenate the gnomAD pass SNVs
# The jobs waits until all bcftools scripts have completed
qsub -b y -cwd -N "concatenate_pass_variants" -e $LOG_DIR -o $LOG_DIR \
    -hold_jid "extract_pass_variants_*" \
    cat "data/scratch/gnomad_v4_*_snps_pass.tsv" > data/interim/gnomad_v4_pass_snvs.tsv