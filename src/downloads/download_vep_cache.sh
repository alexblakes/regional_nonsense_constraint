#!/bin/bash --login

# Download VEP 105 cache.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

URL="https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_vep_105_GRCh38.tar.gz"
DIR="/mnt/bmh01-rds/Ellingford_gene/.vep"
NAME="homo_sapiens_vep_105_GRCh38.tar.gz"
FILEPATH=${DIR}/${NAME}

# Download and extract
wget -P ${DIR} ${URL}
tar xvf ${FILEPATH} --directory ${DIR}