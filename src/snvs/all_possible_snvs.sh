#!/usr/bin/env bash
set -euo pipefail

# Get all possible coding SNVs

BED="data/interim/gencode_v39_canonical_cds.bed"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FILE_OUT="data/interim/cds_all_possible_snvs.vcf.gz"

bedtools makewindows -w 1 -b $BED \
| bedtools getfasta -bedOut -fi $FASTA -bed stdin \
| cut -f 1,3,4 \
| grep -E [ATCG] \
| sort -k1,1V -k2,2n -k3,3 --buffer-size=2G --parallel 8 -u \
| awk \
    -v bases="A,T,C,G" \
    '
        BEGIN {split(bases, b, ",")};
        { for (i in b) { base=b[i]; if($3 != base) { $4 = base; $1 = $1; print $0 } } }
    ' \
| bcftools convert -c CHROM,POS,REF,ALT -f $FASTA --tsv2vcf - -o $FILE_OUT -W=tbi
