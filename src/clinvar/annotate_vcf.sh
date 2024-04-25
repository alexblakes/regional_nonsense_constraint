#!/usr/bin/env bash
# Annotate ClinVar VCF with ACMG and review tags.

VCF="data/interim/clinvar_variants_selected.vcf"
ANNOTS="data/interim/clinvar_variants_selected.tsv"
HEADER="data/manual/header_lines_clinvar_vcf.txt"
FILE_OUT="data/interim/clinvar_variants_annotated.vcf"

# Index annotation file
bgzip --keep --force $ANNOTS
tabix --sequence 1 --begin 2 --end 2 --force "${ANNOTS}.gz"

# Write header file
echo '##INFO=<ID=ACMG,Number=1,Type=String,Description="ACMG classification">' > $HEADER
echo '##INFO=<ID=REVIEW,Number=1,Type=String,Description="ClinVar review status">' >> $HEADER
echo '##INFO=<ID=CV_SYMBOL,Number=1,Type=String,Description="HGNC gene symbol from ClinVar">' >> $HEADER

# Annotate the VCF with ACMG and review tags
bcftools annotate \
    -a "${ANNOTS}.gz" \
    -c CHROM,POS,REF,ALT,CV_SYMBOL,ACMG,REVIEW \
    -h $HEADER \
    $VCF \
> $FILE_OUT

# Clean up
rm "${ANNOTS}.gz" "${ANNOTS}.gz.tbi"