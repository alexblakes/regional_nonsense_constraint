#!/usr/bin/env bash
set -euo pipefail

# Liftover the pext annotations from hg19 to hg38

FILE_IN=data/interim/pext_37.bed
CHAIN=data/raw/hg19ToHg38.over.chain
FILE_OUT=data/interim/pext_38.bed
LOG=data/logs/pext_liftover_unmapped.txt

liftOver $FILE_IN $CHAIN $FILE_OUT $LOG