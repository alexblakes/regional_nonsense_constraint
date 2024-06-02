#!/usr/bin/env bash
set -euo pipefail

# Annotate SNVs

HEADER="data/manual/header_lines_vep_consequence.txt"
POSSIBLE_SNVS="data/interim/cds_all_possible_snvs_vep.tsv.gz"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Write VCF header for VEP consequence
echo '##INFO=<ID=Consequence,Number=1,Type=String,Description="VEP consequence">' > $HEADER

# Annotate SNVs
bcftools convert -c CHROM,POS,REF,ALT,ID,Consequence -f $FASTA --tsv2vcf $POSSIBLE_SNVS \
| bcftools annotate -a $POSSIBLE_SNVS -c CHROM,POS,REF,ALT,~ID,Consequence -h $HEADER \
| grep -v "^##" | less -S