#!/bin/bash --login

# Download reference genome FASTA file.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# FASTA file
wget -c -P data/raw/ https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-tx-annotation/gnomad_browser/all.baselevel.021620.tsv.bgz
gunzip -S .bgz data/raw/all.baselevel.021620.tsv.bgz
