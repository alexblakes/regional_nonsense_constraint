#!/bin/bash --login
set -euo pipefail

# Submit logging script

LOG_DIR="data/logs/csf"
SCRIPT="src/snvs/logging.sh"

qsub \
    -hold_jid "vep_wait_and_concat" \
    -cwd \
    -N "snvs_logging" \
    -e "${LOG_DIR}/\$JOB_NAME.err" \
    -o "${LOG_DIR}/\$JOB_NAME.out" \
    -b y \
    bash $SCRIPT
