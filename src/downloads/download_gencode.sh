#!/bin/bash --login

# Download GENCODE v39 GTF.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz | \
gunzip -c > data/raw/gencode.v39.annotation.gtf
