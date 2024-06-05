#!/bin/bash --login
source activate vep
set -euo pipefail

# Download gnomAD v4.1 coverage data.
# This is apparently identical to the v4.0 coverage data.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# Constants
DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v4.1/exomes/coverage"
URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/coverage/exomes/gnomad.exomes.v4.0.coverage.summary.tsv.bgz"
NAME=$(basename $URL)
OUTFILE=${DIR}/${NAME}

# Make target directory
mkdir -p ${DIR}

# Download and unzip
wget -O - $URL | tee $OUTFILE | md5sum