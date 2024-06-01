#!/usr/bin/env bash

source activate vep # Must precede set command on cluster
set -euo pipefail

# Concatenate VEP-annotated SNVs

TMP="data/tmp"
VCF_PATHS="${TMP}/vcfs.txt"
FILE_OUT="data/interim/cds_all_possible_snvs_vep.vcf.gz"

rm -f $FILE_OUT "${FILE_OUT}.tbi"

find $TMP -type f -name chr*.vcf | sort -k1,1V > $VCF_PATHS
bcftools concat --file-list $VCF_PATHS -W=tbi -o $FILE_OUT

# Clean up 
rm -r $TMP