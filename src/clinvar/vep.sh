#!/usr/bin/env bash

# Annotate ClinVar variants with VEP
# Run in the vep conda environment

FILE_IN="data/interim/clinvar_variants_annotated.vcf.gz"
FILE_OUT="data/interim/clinvar_variants_vep.tsv.gz"
NMD="data/final/nmd_annotations_simple.tsv.gz"
HEADER="data/manual/header_lines_nmd_annotation.txt"

vep \
    --input_file $FILE_IN \
    --species homo_sapiens \
    --assembly GRCh38 \
    --offline \
    --cache \
    --dir_cache /mnt/bmh01-rds/Ellingford_gene/.vep \
    --fork 8 \
    --buffer_size 10000 \
    --force_overwrite \
    --format vcf \
    --vcf \
    --no_stats \
    --minimal \
    --symbol \
    --canonical \
    --fields "Consequence,Feature,SYMBOL,CANONICAL" \
    --output_file - \
| bcftools +split-vep --duplicate --columns - \
| bcftools view -i 'CANONICAL="YES"' \
| bcftools view -i 'Feature=clinvar_enst' \
| bcftools view -i 'Consequence ~ "synonymous\|missense\|stop_gained\|frameshift"' \
| bcftools annotate --set-id %Feature \
| bcftools annotate \
    --remove INFO/CSQ,INFO/CANONICAL,INFO/clinvar_enst,INFO/Feature \
| bcftools annotate -a $NMD -c CHROM,POS,~ID,REGION -h $HEADER \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%ID\t%Consequence\t%REGION\t%ACMG\t%REVIEW' \
| awk -F "\t" -v OFS="\t" \
    '{ \
        if($7 ~ /synonymous/) $7 = "synonymous_variant"; \
        if($7 ~ /missense/) $7 = "missense_variant"; \
        if($7 ~ /stop_gained/) $7 = "stop_gained"; \
        if($7 ~ /frameshift/) $7 = "frameshift_variant"; \
        $1=$1; \
        print $0 \
    }' \
| bgzip > $FILE_OUT