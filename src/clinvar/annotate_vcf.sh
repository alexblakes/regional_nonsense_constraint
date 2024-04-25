#!/usr/bin/env bash
# Annotate ClinVar VCF with ACMG and review tags.

VCF="data/interim/clinvar_variants_selected.vcf"
ANNOTS="data/interim/clinvar_variants_selected.tsv"
HEADER="data/manual/clinvar_vcf_header_lines.txt"
FILE_OUT="data/interim/clinvar_variants_annotated.vcf"

# Index annotation file
bgzip --keep --force $ANNOTS
tabix --sequence 1 --begin 2 --end 2 --force "${ANNOTS}.gz"

# Annotate the VCF with ACMG and review tags
bcftools annotate \
    -a "${ANNOTS}.gz" \
    -c CHROM,POS,REF,ALT,CV_SYMBOL,ACMG,REVIEW \
    -h $HEADER \
    $VCF \
> $FILE_OUT

# Clean up
rm "${ANNOTS}.gz" "${ANNOTS}.gz.tbi"