#!/usr/bin/env bash
set -euo pipefail

GENE=$1
CHR=$2
START=$3
END=$4

CLINVAR="data/interim/clinvar_variants_annotated.vcf"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FILE_TMP="data/interim/clinvar_vep_protein_paint_${GENE}.vcf"
FILE_OUT_VUS="data/final/protein_paint/${GENE}_clinvar_vus.vcf.gz"
FILE_OUT_PLP="data/final/protein_paint/${GENE}_clinvar_plp.vcf.gz"

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
    --fields "Consequence,Feature,SYMBOL,CANONICAL,HGVSp,Protein_position,Amino_acids,HGVSc,Existing_variation" \
    --output_file STDOUT \
| filter_vep \
    --filter "SYMBOL is $GENE" \
    --filter "CANONICAL is YES" \
    --filter "Consequence in stop_gained,frameshift_variant" \
    --only_matched \
> $FILE_TMP

bcftools view -i 'ACMG = "VUS"' -Oz --write-index=tbi -o $FILE_OUT_VUS $FILE_TMP
bcftools view -i 'ACMG ~ "LP\|P"' -Oz --write-index=tbi -o $FILE_OUT_PLP $FILE_TMP
