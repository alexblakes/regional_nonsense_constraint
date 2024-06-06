#!/bin/bash --login

# Download *indexed* VEP 105 cache.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

URL="https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz"
DIR="/mnt/bmh01-rds/Ellingford_gene/.vep"
NAME=$(basename $URL)
FILEPATH="${DIR}/${NAME}"

# Download and extract
wget -O - ${URL} | tee $FILEPATH | sum -
tar xvf ${FILEPATH} --directory ${DIR}