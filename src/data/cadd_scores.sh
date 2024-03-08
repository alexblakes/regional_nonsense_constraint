#!/usr/bin/env bash

# Get the lowest CADD score for each genomic position
tabix data/raw/whole_genome_SNVs.tsv.gz \
    -R data/interim/gencode_v39_canonical_cds_no_chr.bed \
    --cache 10000 \
| cut -f 1,2,6 \
| sort -k1,1V -k2,2n -k3,3n -S128G --parallel 8 \
| sort -u -k1,1V -k2,2n -s -S128G --parallel 8 \
| sed 's/^/chr/' \
> data/interim/cadd_all_sites_min_scores.tsv