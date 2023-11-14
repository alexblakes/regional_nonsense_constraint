#!/bin/bash --login

# Concatenate the gnomAD SNVs
# The job waits until all bcftools scripts have completed

#$ -cwd
#$ -hold_jid "extract_pass_variants_*"
#$ -N "concatenate_pass_variants"
#$ -e data/logs/csf/
#$ -o data/logs/csf/

cat data/scratch/gnomad_v4_*_snps_pass.tsv > data/interim/gnomad_v4_pass_snvs.tsv