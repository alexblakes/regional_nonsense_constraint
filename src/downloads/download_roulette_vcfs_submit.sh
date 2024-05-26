#!/bin/bash --login

set -euo pipefail

# Submission script to download Roulette data.
# Run on the login node
# Run from ukb_constraint directory

CHRS=({1..22}) # PRODUCTION

for CHR in ${CHRS[@]};
    do
        VCF_URL="http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/${CHR}_rate_v5.2_TFBS_correction_all.vcf.bgz"
        qsub \
            -b y \
            -o data/logs/csf/ \
            -e data/logs/csf/ \
            -cwd \
            /mnt/bmh01-rds/Ellingford_gene/ukb_constraint/src/downloads/download_roulette_vcfs.sh ${VCF_URL}
    done