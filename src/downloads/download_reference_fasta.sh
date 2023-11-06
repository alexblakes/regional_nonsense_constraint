#!/bin/bash --login

# Download reference genome FASTA file.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# FASTA file
wget -c -P data/raw/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# FASTA index
wget -c -P data/raw/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai