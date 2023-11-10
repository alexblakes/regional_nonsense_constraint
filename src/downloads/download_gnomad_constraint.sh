#!/bin/bash --login

# Download gnomAD v4.0 constraint data.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# Constants
DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v4.0/constraint"
URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/v4.0/constraint/gnomad.v4.0.constraint_metrics.tsv"
README_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/v4.0/constraint/README.txt"

# Download
wget -c -P ${DIR} ${URL}
wget -c -P ${DIR} ${README_URL}