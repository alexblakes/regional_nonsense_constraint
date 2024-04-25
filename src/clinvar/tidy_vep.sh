#!/usr/bin/env bash
# Filter VEP-annotated ClinVar variants.

FILE_IN="data/interim/clinvar_variants_vep.vcf"
FILE_OUT="data/interim/clinvar_variants_vep_tidy.tsv"
NMD="data/final/nmd_annotations_simple.tsv.gz"
HEADER="data/manual/header_lines_nmd_annotation.txt"

bcftools +split-vep $FILE_IN --duplicate --columns - \
| bcftools view -i 'CANONICAL="YES"' \
| bcftools view -i 'CV_SYMBOL=SYMBOL' \
| bcftools view -i 'Consequence ~ "synonymous\|missense\|stop_gained\|frameshift"' \
| bcftools annotate \
    --set-id '%Feature' \
| bcftools annotate \
    --remove INFO/CSQ,INFO/CANONICAL,INFO/CV_SYMBOL,INFO/SYMBOL,INFO/REVIEW,INFO/Feature \
| bcftools annotate -a $NMD -c CHROM,POS,~ID,REGION -h $HEADER \
| bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%Consequence\t%REGION\t%ACMG' \
| awk -F "\t" -v OFS="\t" \
    '{ \
        if($6 ~ /synonymous/) $6 = "synonymous_variant"; \
        if($6 ~ /missense/) $6 = "missense_variant"; \
        if($6 ~ /stop_gained/) $6 = "stop_gained"; \
        if($6 ~ /frameshift/) $6 = "frameshift_variant"; \
        $1=$1; \
        print $0 \
    }' \
> $FILE_OUT
