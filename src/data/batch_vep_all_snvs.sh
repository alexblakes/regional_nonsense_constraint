#!/bin/bash --login

# Submits jobs to annotate SNVs with VEP.

# Constants
LOG_DIR="data/logs/csf/"

# Submit one job per chromosome
for N in {0..29}; do
  qsub -b y -cwd -N "vep_snvs_${N}" -e $LOG_DIR -o $LOG_DIR \
  src/data/vep_all_snvs.sh $N
done