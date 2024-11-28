#!/usr/bin/env bash
# Group genes into paralog families
set -euo pipefail
source src/utils.sh

FILE_LOG="data/logs/src.paralogues.get_paralog_families.log"
FILE_IN="data/interim/ensembl_paralogs_filtered.tsv"
DIR_TMP=$(mkdir -p data)

trap 'rm -rf -- "${DIR_TMP}"' EXIT

< $FILE_IN less -S