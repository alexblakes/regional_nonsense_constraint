#!/usr/bin/env/ bash

# This script extracts all CDS SNVs passing gnomAD filters using bcftools

# A function to run bcftools
function bcf {
    bcftools view \
        -v snps \
        -f PASS \
        --min-ac 1 \
        -R "data/interim/gencode_v39_canonical_cds.bed" \
        -O u \
        "${1}" |\
    bcftools head \
        -n 10 \
        -O u |\
    bcftools query \
        -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t.\t.\t%AC\t%AN\n' # \
        # -o "$(basename "${1}" .vcf.gz)_snps_pass.tsv"
    echo "Done with ${1} at $(date +"%T")"
}

bcf $1

# # Parallel processing

# ## Export function for use by parallel
# export -f bcf

# ## Install parallel
# apt-get update
# apt-get install parallel -y

# ## Run jobs in parallel
# find vcfs/* | grep -v .tbi | parallel bcf

# # Print finish notice
# echo "Finished with $(find split_paths_*) at $(date +"%T")"