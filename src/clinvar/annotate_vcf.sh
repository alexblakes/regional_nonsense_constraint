#!/usr/bin/env bash
# Annotate ClinVar VCF with ACMG and review tags.

VCF="data/interim/clinvar_variants_selected.vcf.gz"
ANNOTS="data/interim/clinvar_variants_selected.tsv"
HEADER="data/manual/header_lines_clinvar_vcf.txt"
FILE_OUT="data/interim/clinvar_variants_annotated.vcf.gz"

# Index annotation file
< $ANNOTS tail -n+2 | sort -k1,1V -k2,2n | bgzip > "${ANNOTS}.gz"
tabix --sequence 1 --begin 2 --end 2 --force "${ANNOTS}.gz"

# Write header file
echo '##INFO=<ID=ACMG,Number=1,Type=String,Description="ACMG classification">' > $HEADER
echo '##INFO=<ID=REVIEW,Number=1,Type=String,Description="ClinVar review status">' >> $HEADER
echo '##INFO=<ID=clinvar_enst,Number=1,Type=String,Description="Ensembl transcript ID from ClinVar">' >> $HEADER

# Annotate the VCF with ACMG and review tags
< $VCF bcftools annotate \
    -a "${ANNOTS}.gz" \
    -c CHROM,POS,REF,ALT,-,ACMG,REVIEW,-,clinvar_enst \
    -h $HEADER \
    -o $FILE_OUT \
    -W=tbi