#!/usr/bin/env bash
# Convert filtered ClinVar variants to VCF.

FILE_IN="data/interim/clinvar_variants_selected.tsv"
FILE_OUT="data/interim/clinvar_variants_selected.vcf.gz"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

tail -n +2 $FILE_IN \
| bcftools convert \
    -c CHROM,POS,REF,ALT \
    -f $FASTA \
    --tsv2vcf - \
| bcftools sort -o $FILE_OUT -W=tbi
