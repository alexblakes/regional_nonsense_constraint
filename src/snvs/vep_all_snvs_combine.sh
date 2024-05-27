#!/usr/bin/env bash

source activate vep # Must precede set command on cluster
set -euo pipefail

# Concatenate VEP-annotated SNVs

TMP="data/tmp"
VCF_PATHS="${TMP}/vcfs.txt"
FILE_OUT="data/interim/cds_all_possible_snvs_vep.vcf.gz"

find $TMP -type f -name chr*.vcf.gz | sort -k1,1V > $VCF_PATHS
bcftools concat --file-list $VCF_PATHS --naive -o $FILE_OUT
tabix -f $FILE_OUT

# Clean up 
rm -r $TMP