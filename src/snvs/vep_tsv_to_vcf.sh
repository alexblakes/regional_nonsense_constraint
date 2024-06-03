#!/usr/bin/env bash
set -euo pipefail

# Convert combined VEP output to VCF

HEADER="data/manual/header_lines_vep_consequence.txt"
VEP_TSV="data/interim/cds_all_possible_snvs_vep.tsv.gz"
VEP_VCF="data/interim/cds_all_possible_snvs_vep.vcf.gz"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Write VCF header for VEP consequence
echo '##INFO=<ID=Consequence,Number=1,Type=String,Description="VEP consequence">' > $HEADER

# Convert VEP output to indexed VCF
bcftools convert \
    -c CHROM,POS,REF,ALT,ID,Consequence \
    -f $FASTA \
    --tsv2vcf $VEP_TSV \
    -W=tbi \
    -o $VEP_VCF
