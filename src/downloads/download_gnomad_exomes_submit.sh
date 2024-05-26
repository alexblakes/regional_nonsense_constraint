#!/bin/bash --login

# Submission script to download gnomAD data.
# Run on the login node
# Run from ukb_constraint directory

CHRS=({1..22} X Y) # PRODUCTION

for CHR in ${CHRS[@]};
    do
        TBI_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr${CHR}.vcf.bgz.tbi"
        VCF_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr${CHR}.vcf.bgz"
        qsub \
            -b y \
            -o data/logs/csf/ \
            -e data/logs/csf/ \
            -cwd \
            /mnt/bmh01-rds/Ellingford_gene/ukb_constraint/src/downloads/download_gnomad_exomes.sh ${TBI_URL} ${VCF_URL}
    done