#!/usr/bin/env bash
set -euo pipefail

GENE=$1
CHR=$2
START=$3
END=$4

CLINVAR="data/interim/clinvar_variants_annotated.vcf"
FILE_OUT="data/interim/clinvar_${GENE}.vcf.gz"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

bgzip -kf $CLINVAR
tabix -f "${CLINVAR}.gz"

bcftools view \
    -r "${CHR}:${START}-${END}" "${CLINVAR}.gz" \
    -O v \
    -o - \
|vep \
    --species homo_sapiens \
    --assembly GRCh38 \
    --fasta $FASTA \
    --offline \
    --cache \
    --dir_cache /mnt/bmh01-rds/Ellingford_gene/.vep \
    --fork 8 \
    --buffer_size 50000 \
    --force_overwrite \
    --format vcf \
    --vcf \
    --no_stats \
    --minimal \
    --symbol \
    --canonical \
    --hgvs \
    --fields "Consequence,Feature,SYMBOL,CANONICAL,HGVSp" \
    --output_file STDOUT \
| filter_vep \
    --filter "SYMBOL is $GENE" \
    --filter "CANONICAL is YES" \
    --filter "Consequence in stop_gained,frameshift_variant" \
    --only_matched
    # -o $FILE_OUT
# > $FILE_OUT
