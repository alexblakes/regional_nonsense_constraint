#!/bin/bash --login
set -euo pipefail

# Submit jobs to annotate SNVs with VEP.

LOG_DIR="data/logs/csf"

qsub \
    -cwd \
    -N "roulette_scaling" \
    -e "${LOG_DIR}/\$JOB_NAME.err" \
    -o "${LOG_DIR}/\$JOB_NAME.out" \
    -b y \
    bash src/roulette/add_scaled_rates.sh