#!/bin/bash

# Extract SNVs passing gnomAD filters using bcftools

# Load conda environment for bcftools
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate
conda activate bio

# A function to run bcftools
function bcf {
    bcftools view \
        -v snps \
        -f PASS \
        --min-ac 1 \
        -R "data/interim/gencode_v39_canonical_cds.bed" \
        -O u \
        "/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v4.0/vcf/gnomad.exomes.v4.0.sites.${1}.vcf.bgz" |\
    bcftools query \
        -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t.\t.\t%AC\t%AN\n' \
        -o "data/scratch/gnomad_v4_${1}_snps_pass.tsv"
    echo "Done with ${1} at $(date +"%T")"
}

# Execute the function
bcf $1