#!/usr/bin/env bash

source activate vep # Must precede the 'set' command on the cluster
set -euo pipefail

# Concatenate VEP-annotated SNVs

TMP="data/tmp"
VCF_PATHS="${TMP}/vcfs.txt"
FILE_OUT="data/interim/cds_all_possible_snvs_vep.tsv.gz"

# Write the filepaths of the intermediate VCFs to a text file
find $TMP -type f -name chr*.vcf | sort -k1,1V > $VCF_PATHS

# Concatenate the VCFs.
# Sanitise the consequences (e.g. drop splicing consequences).
bcftools concat --file-list $VCF_PATHS \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%Consequence' \
| awk -v OFS="\t" \
    '
        $6 ~ /synonymous/ {$6="synonymous_variant"}
        $6 ~ /missense/ {$6="missense_variant"}
        $6 ~ /stop_gained/ {$6="stop_gained"}
        {$1=$1; print $0}
    ' \
| bgzip > $FILE_OUT
tabix -f -s1 -b2 -e2 $FILE_OUT

# Clean up 
wait
rm -r $TMP