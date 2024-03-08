#!/bin/bash --login

# Download CADD v1.7 data.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# Constants
DIR="data/raw/"
URL_TBI="https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi"
URL_GZ="https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"

# Download
wget -c -P ${DIR} --no-check-certificate ${URL_TBI}
wget -c -P ${DIR} --no-check-certificate ${URL_GZ}