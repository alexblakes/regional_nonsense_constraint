#!/usr/bin/env bash
# Filter VEP-annotated ClinVar variants.

FILE_IN="data/interim/clinvar_variants_vep.vcf"
FILE_OUT="data/interim/clinvar_variants_vep_clean.vcf"
NMD="data/final/nmd_annotations_simple.tsv.gz"
HEADER="data/manual/header_lines_nmd_annotation.txt"

bcftools +split-vep $FILE_IN --duplicate --columns - \
| bcftools view -i 'CANONICAL="YES"' \
| bcftools view -i 'CV_SYMBOL=SYMBOL' \
| bcftools view -i 'Consequence ~ "synonymous\|missense\|stop_gained"' \
| bcftools annotate \
    --set-id '%Feature' \
| bcftools annotate \
    --remove INFO/CSQ,INFO/CANONICAL,INFO/CV_SYMBOL,INFO/SYMBOL,INFO/REVIEW,INFO/Feature \
| bcftools annotate -a $NMD -c CHROM,POS,~ID,REGION -h $HEADER \
| grep -v "^#" | less -S
# -i INFO/CV_SYMBOL=INFO/SYMBOL \
# --columns Consequence,Feature,SYMBOL,CANONICAL \