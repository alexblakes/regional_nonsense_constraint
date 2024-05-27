#!/bin/bash --login
# Submit jobs to annotate SNVs with VEP.

set -euo pipefail

LOG_DIR="data/logs/csf"
TMP="data/tmp"
VCF_PATHS="${TMP}/vcfs.txt"
FILE_OUT="data/interim/cds_all_possible_snvs_vep.vcf.gz"

mkdir -p $TMP

# Submit one job per chromosome
for c in chr{1..22}; do
  qsub \
    -cwd \
    -N "vep_snvs_${c}" \
    -e "${LOG_DIR}/\$JOB_NAME.err" \
    -o "${LOG_DIR}/\$JOB_NAME.out" \
    -b y \
  bash src/snvs/vep_all_snvs.sh $c
done

# Wait for jobs to complete, then combine VCFs
qsub \
    -hold_jid "vep_snvs_chr*" \
    -cwd \
    -N "vep_wait_and_concat" \
    -e "${LOG_DIR}/\$JOB_NAME.err" \
    -o "${LOG_DIR}/\$JOB_NAME.out" \
    -b y \
    bash src/snvs/vep_all_snvs_combine.sh