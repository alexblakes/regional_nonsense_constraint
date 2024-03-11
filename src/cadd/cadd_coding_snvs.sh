#! /usr/bin/env/ bash

# Subset CADD data to coding SNVs only

tabix data/raw/whole_genome_SNVs.tsv.gz \
    -R data/interim/gencode_v39_canonical_cds_no_chr.bed \
    > data/interim/cadd_scores_coding.tsv