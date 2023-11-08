#!/bin/bash --login

# Submission script to download gnomAD data.
# Run from login node

CHRS=({1..22} X Y)
# CHRS=({1..2}) # TESTING

for CHR in ${CHRS[@]};
    do
        TBI_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr"${CHR}".vcf.bgz.tbi"
        VCF_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr"${CHR}".vcf.bgz"
        qsub -b y ~/scratch/ukb_constraint/src/downloads/download_gnomad.sh ${TBI_URL} ${VCF_URL}
        # bash ~/scratch/ukb_constraint/src/downloads/download_gnomad.sh ${TBI_URL} ${VCF_URL}
    done