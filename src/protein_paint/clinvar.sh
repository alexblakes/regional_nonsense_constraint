#!/usr/bin/env bash
set -euo pipefail

GENE=$1
CHR=$2
START=$3
END=$4

CLINVAR="data/interim/clinvar_variants_annotated.vcf.gz"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FILE_TMP="data/interim/clinvar_vep_protein_paint_${GENE}.vcf"
FILE_OUT_VUS="data/final/protein_paint/${GENE}_clinvar_vus.txt"
FILE_OUT_PLP="data/final/protein_paint/${GENE}_clinvar_plp.txt"

bcftools view \
    -r "${CHR}:${START}-${END}" "${CLINVAR}" \
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
    --canonical \
    --hgvs \
    --fields "Consequence,Feature,SYMBOL,CANONICAL,HGVSp,HGVSc,Protein_position" \
    --output_file STDOUT \
| filter_vep \
    --filter "SYMBOL is $GENE" \
    --filter "CANONICAL is YES" \
    --filter "Consequence in stop_gained,frameshift_variant" \
    --only_matched \
| bcftools +split-vep \
    --columns - \
    --format '%HGVSc\t%Protein_position\t%Consequence\t%ACMG' \
| awk -v OFS="\t" \
    ' \
        { \
            sub(/.*:c\./, "", $1) \
            sub(/\-[0-9]+/, "", $2) \
            sub("frameshift_variant", "F", $3) \
            sub("stop_gained", "N", $3); \
            $1=$1; print $0 \
        } \
    ' \
> $FILE_TMP

awk -v OFS="\t" '$NF == "VUS" {NF-=1; print $0}' $FILE_TMP > $FILE_OUT_VUS
awk -v OFS="\t" '$NF ~ /P|LP/ {NF-=1; print $0}' $FILE_TMP > $FILE_OUT_PLP

# Clean up
rm $FILE_TMP
