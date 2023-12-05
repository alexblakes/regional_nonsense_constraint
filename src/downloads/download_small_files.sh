#!/bin/bash --login

# Download small files.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# Cardiac G2P
wget -O - https://www.ebi.ac.uk/gene2phenotype/downloads/CardiacG2P.csv.gz | \
gunzip -c > data/raw/CardiacG2P.csv

# Skeletal G2P
wget -O - https://www.ebi.ac.uk/gene2phenotype/downloads/SkeletalG2P.csv.gz | \
gunzip -c > data/raw/SkeletalG2P.csv

# DDG2P
wget -O - https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz | \
gunzip -c > data/raw/DDG2P.csv

# Eye G2P
wget -O - https://www.ebi.ac.uk/gene2phenotype/downloads/EyeG2P.csv.gz | \
gunzip -c > data/raw/EyeG2P.csv

# SkinG2P
wget -O - https://www.ebi.ac.uk/gene2phenotype/downloads/SkinG2P.csv.gz | \
gunzip -c > data/raw/SkinG2P.csv

# Liftover chain
wget -O - http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz | \
gunzip -c > data/raw/hg19ToHg38.over.chain

# ClinVar variant summary
wget -O - https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz | \
gunzip -c > data/raw/variant_summary.txt

# OMIM genemap2
wget -P data/raw/ -c --no-check-certificate https://data.omim.org/downloads/dW5qBh3GSkCt5K3BqUQL1w/genemap2.txt

# gnomAD v4.0 constraint
wget -P data/raw/ https://storage.googleapis.com/gcp-public-data--gnomad/release/v4.0/constraint/gnomad.v4.0.constraint_metrics.tsv

# gnomAD non-coding mutation rates
wget -P data/raw/ https://storage.googleapis.com/gnomad-nc-constraint-v31-paper/mutation_rate_by_context_methyl.txt

# SNV-level Alpha Missense scores
wget -O - https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz | \
gunzip -c > data/raw/AlphaMissense_hg38.tsv
