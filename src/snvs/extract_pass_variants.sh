#!/usr/bin/env bash
set -euo pipefail

# Extract SNVs on autosomes passing gnomAD filters using bcftools

export REGIONS="data/interim/gencode_v39_canonical_cds.bed"
export GNOMAD_DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v4.1/exomes/vcf"
export TMP_DIR="data/tmp"
VCF_LIST="${TMP_DIR}/vcf_list.txt"
FILE_OUT="data/interim/gnomad_v4.1_pass_snvs.vcf.gz"

mkdir -p $TMP_DIR

extract_pass_snvs() {

    local VCF="gnomad.exomes.v4.1.sites.${1}.vcf.bgz"
    local FILE_IN="${GNOMAD_DIR}/${VCF}"
    local BASENAME=$(basename -s .vcf.bgz $VCF)
    local FILE_OUT="${TMP_DIR}/${BASENAME}.snvs.pass.vcf.gz"

    bcftools view \
        --type snps \
        --apply-filters PASS \
        --min-ac 1 \
        -R $REGIONS \
        -T $REGIONS \
        -Ou \
        $FILE_IN \
    | bcftools annotate \
        -x "^INFO/AC,INFO/AN,INFO/AF" \
        --set-id "." \
        -o $FILE_OUT

    echo "Done with ${1} at $(date +"%T")"
}

# Run in each chromosome in parallel
export -f extract_pass_snvs
echo chr{1..22} | xargs -n1 | parallel extract_pass_snvs {}

# Concatenate the outputs
find $TMP_DIR -type f -name *.vcf.gz \
| sort -k1,1V \
> $VCF_LIST

bcftools concat -f $VCF_LIST --naive > $FILE_OUT
tabix $FILE_OUT

# Clean up
rm -r $TMP_DIR
