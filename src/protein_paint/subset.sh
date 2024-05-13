#!/usr/bin/env bash
set -euo pipefail

GENE=$1
CHR=$2
START=$3
END=$4

CLINVAR="data/interim/clinvar_variants_annotated.vcf"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FILE_TMP="data/interim/clinvar_vep_protein_paint_${GENE}.vcf"
FILE_OUT_VUS="data/final/protein_paint/${GENE}_clinvar_vus.txt"
FILE_OUT_PLP="data/final/protein_paint/${GENE}_clinvar_plp.txt"
HEADER="gene\trefseq\tchromosome\tstart\taachange\tclass"

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
    --merged \
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
    --fields "SOURCE,Consequence,Feature,SYMBOL,CANONICAL,HGVSp,HGVSc" \
    --output_file STDOUT \
| filter_vep \
    --filter "SOURCE is RefSeq" \
    --filter "SYMBOL is $GENE" \
    --filter "CANONICAL is YES" \
    --filter "Consequence in stop_gained,frameshift_variant" \
    --only_matched \
| bcftools +split-vep \
    --columns - \
    --format '%CHROM\t%POS\t%REF\t%ALT\t%ACMG\t%Consequence\t%Feature\t%SYMBOL\t%HGVSp' \
| awk -v OFS="\t" \
    ' \
        { \
            sub(/\.[0-9]+/, "", $7); \
            sub("frameshift_variant", "frameshift", $6); \
            sub("stop_gained", "nonsense", $6); \
            print($8, $7, $1, $2, $NF, $6, $5) \
        } \
    ' \
> $FILE_TMP

awk -v OFS="\t" -v H=$HEADER 'BEGIN {print H}; $NF == "VUS" {NF-=1; print $0}' $FILE_TMP > $FILE_OUT_VUS

# Clean up
rm $FILE_TMP