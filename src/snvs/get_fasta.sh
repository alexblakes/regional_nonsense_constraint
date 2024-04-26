#!/usr/bin/env bash

# Get the FASTA sequences for canonical CDS features

# The coordinates of each region are extended by 1nt in either direction. This is so 
# that the sequence context for the most extreme 5' and 3' CDS positions can still be
# obtained.

FILE_IN="data/interim/gencode_v39_canonical_cds.bed"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FILE_OUT="data/interim/gencode_v39_canonical_cds_seq.tsv"

bedtools slop \
    -i $FILE_IN \
    -g "${FASTA}.fai" \
    -b 1 |\
bedtools getfasta \
    -fo $FILE_OUT \
    -tab \
    -fi $FASTA \
    -bed stdin

