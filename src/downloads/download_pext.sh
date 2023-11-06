#!/bin/bash --login

# Download reference genome FASTA file.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# FASTA file
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-tx-annotation/gnomad_browser/all.baselevel.021620.tsv.bgz | \
gunzip -c -S .bgz > data/raw/all.baselevel.021620.tsv
