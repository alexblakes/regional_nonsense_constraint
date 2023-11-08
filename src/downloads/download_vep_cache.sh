#!/bin/bash --login

# Download vep 105 cache.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

DIR="~/mnt/bmh01-rds/Ellingford_gene/.vep"

# FASTA file
wget -P ${DIR} https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_vep_105_GRCh38.tar.gz
tar xzf ""${DIR}"/homo_sapiens_vep_105_GRCh38.tar.gz" --directory ${DIR}