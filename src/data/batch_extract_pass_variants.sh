#!/bin/bash --login

# Submits jobs which extract gnomAD variants to the CSF.
CHRS=( {1..22} X Y )
CHRS=( "${CHRS[@]/#/chr}" )

LOG_DIR="data/logs/csf/"

for CHR in ${CHRS[@]}; do
  qsub -b y -cwd -N "extract_pass_variants_${CHR}" -e $LOG_DIR -o $LOG_DIR \
  src/data/extract_pass_variants.sh $CHR
done