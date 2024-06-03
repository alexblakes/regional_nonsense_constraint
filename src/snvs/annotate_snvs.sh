#!/usr/bin/env bash
set -euo pipefail

# Annotate SNVs

HEADER_CSQ="data/manual/header_lines_vep_consequence.txt"
VEP_TSV="data/interim/cds_all_possible_snvs_vep.tsv.gz"
VEP_VCF="data/interim/cds_all_possible_snvs_vep.vcf.gz"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
NMD="data/final/nmd_annotations_simple.tsv.gz"
HEADER_NMD="data/manual/header_lines_nmd_annotation.txt"
GNOMAD_GZ="data/interim/gnomad_v4.1_pass_snvs.vcf.gz"

# The -a file must be in tsv format in order for bcftools annotate to work in a stream.

# Annotate SNVs
bcftools annotate -a $GNOMAD_GZ -c AC,AN,AF $VEP_VCF -Oz \
| bcftools annotate -a $VEP_VCF -c ~ID,Consequence \
| grep -v "^##" | less -S
# | bcftools annotate -a $NMD -c CHROM,POS,~ID,REGION -h $HEADER_NMD -Ou \