#!/usr/bin/env bash
set -euo pipefail

FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FILE_IN="data/interim/cds_all_possible_snvs.vcf"
FILE_OUT="data/interim/cds_all_possible_snvs.vcf.gz"

bcftools convert -c CHROM,POS,ID,REF,ALT,-,-,- -f $FASTA --tsv2vcf $FILE_IN | less -S