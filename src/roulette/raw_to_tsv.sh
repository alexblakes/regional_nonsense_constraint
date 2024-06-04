#!/usr/bin/env bash
set -euo pipefail

# Convert raw Roulette scores to TSV for annotation
## Where an adjusted Roulette score is available, use it instead of the raw score.

FILE_IN="data/interim/roulette/raw_cds_combined.vcf.gz"
FILE_OUT="data/interim/roulette/raw_cds_combined.tsv.gz"

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%PN\t%MR\t%AR\t%MG\t%MC' \
    $FILE_IN \
| sed 's/^/chr/' \
| awk -v OFS="\t" '
    {$(NF+1) = $7};
    $8 != "." {$NF = $8};
    {$1=$1; print $0}
' \
| sort -k1,1V -k2,2n -k4,4 -u --parallel 8 --buffer-size 20% \
| bgzip -f > $FILE_OUT

tabix -s1 -b2 -e2 $FILE_OUT