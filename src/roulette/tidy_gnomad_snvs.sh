#!/usr/bin/env bash
set -euo pipefail

# Tidy the observed SNVs in gnomAD for Roulette

FILE_IN="data/interim/gnomad_v4.1_pass_snvs.vcf.gz"
FILE_OUT="data/interim/gnomad_v4.1_pass_snvs_roulette.tsv"

( printf "CHROM\tPOS\tREF\tALT\n"; \
bcftools view $FILE_IN \
| grep -v "^#" \
| cut -f1,2,4,5 \
| sed 's/chr//') \
> $FILE_OUT