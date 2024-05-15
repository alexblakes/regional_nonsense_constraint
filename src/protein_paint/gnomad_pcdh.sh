#!/usr/bin/env bash
set -euo pipefail

# Extract variants from gnomAD for the protocadherin gene cluster

GENE=$1
CHR=$2
START=$3
END=$4
CSQ=$5

VCF_DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v4.0/vcf/"
VCF_FILE="gnomad.exomes.v4.0.sites.${CHR}.vcf.bgz"
VCF="${VCF_DIR}${VCF_FILE}"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FILE_OUT="data/final/protein_paint/${GENE}_gnomad_${CSQ}.txt"

bcftools view \
    -r "${CHR}:${START}-${END}" $VCF \
    -i 'FILTER="PASS" & AC>=1' \
    -O v \
    -o - \
|vep \
    --species homo_sapiens \
    --assembly GRCh38 \
    --fasta $FASTA \
    --offline \
    --cache \
    --dir_cache /mnt/bmh01-rds/Ellingford_gene/.vep \
    --fork 4 \
    --force_overwrite \
    --format vcf \
    --vcf \
    --no_stats \
    --minimal \
    --symbol \
    --pick \
    --hgvs \
    --fields "Consequence,Feature,SYMBOL,CANONICAL,HGVSc,Protein_position" \
    --output_file STDOUT \
| filter_vep \
    --filter "Consequence is ${CSQ}" \
    --only_matched \
| bcftools +split-vep \
    --columns - \
    --format '%HGVSc\t%POS\t%Consequence\t%AC' \
| awk -v OFS="\t" \
    ' \
        { \
            sub(/.*:c\./, "", $1); \
            sub(/-[0-9]+/, "", $2); \
            sub("stop_gained", "N", $3); \
            sub("frameshift_variant", "F", $3); \
            $1=$1; \
            N = 1; \
            while (N <= $4) { \
                print $1, $2, $3; \
                N += 1 ; \
            } \
        } \
    ' \
> $FILE_OUT