#!/bin/bash
set -euo pipefail

# Script to download Roulette VCF for one chromosome.

_DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs"
_FILE=$(basename $1)

wget -c -P ${_DIR} $1

md5sum "${_DIR}/${_FILE}"
