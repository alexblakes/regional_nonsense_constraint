#!/usr/bin/env bash
set -euo pipefail

# Convert the gnomAD PASS SNVs to TSV

FILE_IN="data/interim/gnomad_v4.1_pass_snvs.vcf.gz"
HEADER="data/manual/header_lines_gnomad.txt"
FILE_OUT="data/interim/gnomad_v4.1_pass_snvs.tsv.gz"

# Write header lines for later annotation
cat << EOF > $HEADER
##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count of PASS SNVs in gnomAD 4.1">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number of PASS SNVs in gnomAD 4.1">
##INFO=<ID=AF,Number=1,Type=String,Description="Allele frequency of PASS SNVs in gnomAD 4.1">
EOF

# Convert to TSV
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF' $FILE_IN \
| bgzip -c > $FILE_OUT && tabix -s1 -b2 -e2 $FILE_OUT