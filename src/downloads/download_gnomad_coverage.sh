#!/bin/bash --login

# Download gnomAD v4.0 coverage data.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# Constants
DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v4.0/coverage"
URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/coverage/exomes/gnomad.exomes.v4.0.coverage.summary.tsv.bgz"
NAME="gnomad.exomes.v4.0.coverage.summary.tsv"
PATH=${DIR}/${NAME}

# Make target directory
mkdir -p ${DIR}

# Download and unzip
wget -O - $URL | gunzip -c -S .bgz > ${PATH}

# MD5 check
md5sum ${PATH} > ${DIR}/coverage.md5