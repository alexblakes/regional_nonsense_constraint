#!/bin/bash --login

# Download reference genome FASTA file.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# Cardiac G2P
wget -P data/raw/ -c https://www.ebi.ac.uk/gene2phenotype/downloads/CardiacG2P.csv.gz
gunzip data/raw/CardiacG2P.csv.gz

# Skeletal G2P
wget -P data/raw/ -c https://www.ebi.ac.uk/gene2phenotype/downloads/SkeletalG2P.csv.gz
gunzip data/raw/SkeletalG2P.csv.gz

# DDG2P
wget -P data/raw/ -c https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz
gunzip data/raw/DDG2P.csv.gz

# Eye G2P
wget -P data/raw/ -c https://www.ebi.ac.uk/gene2phenotype/downloads/EyeG2P.csv.gz
gunzip data/raw/EyeG2P.csv.gz

# SkinG2P
wget -P data/raw/ -c https://www.ebi.ac.uk/gene2phenotype/downloads/SkinG2P.csv.gz
gunzip data/raw/SkinG2P.csv.gz

# Liftover chain
wget -P data/raw/ -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip data/raw/hg19ToHg38.over.chain.gz

# ClinVar variant summary
wget -P data/raw/ -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip data/raw/variant_summary.txt.gz

# OMIM genemap2
wget -P data/raw/ -c --no-check-certificate https://data.omim.org/downloads/dW5qBh3GSkCt5K3BqUQL1w/gendata/raw/emap2.txt
