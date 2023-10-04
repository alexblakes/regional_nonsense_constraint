#!/usr/bin/env/ bash

# This script extracts all CDS SNVs failing gnomAD filters using bcftools

# Path to example WGS pVCF
vcf="/mnt/project/Bulk/Exome sequences_Alternative exome processing/Exome variant call files (gnomAD) (VCFs)/ukb24068_c10_b1002_v1.vcf.gz"

# A function to run bcftools
function bcf {
    bcftools view \
        -v snps \
        -e 'FILTER="PASS"' \
        -R "./gencode_v39_canonical_cds_chr.bed" \
        -O u \
        "${1}" |\
    bcftools query \
        -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t%FILTER\t.\t%AC\t%AN\n' \
        -o "$(basename "${1}" .vcf.gz)_snps_fail.tsv"
    echo "Done with ${1} at $(date +"%T")"
}

# Download VCFs and indexes
rm -rf vcfs/
mkdir vcfs/
xargs -a split_paths_* -d "\n" dx download -o vcfs/

# Parallel processing

## Export function for use by parallel
export -f bcf

## Install parallel
apt-get update
apt-get install parallel -y

## Run jobs in parallel
find vcfs/* | grep -v .tbi | parallel bcf

# Tidy working directory
rm -r vcfs/

# Print finish notice
echo "Finished with $(find split_paths_*) at $(date +"%T")"