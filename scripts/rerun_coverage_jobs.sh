#!/usr/bin/env bash

# This script reruns coverage jobs which have failed

# Change these variables for each batch of jobs
JOB_NAME="*production_230408*"
TAG_01="original_01"
RERUN="rerun_01"
instance_type="mem3_ssd1_v2_x4"
instances=16
cost_limit=1.50

# No need to change these variables
JOB_ID_FILE="failed_jobs.txt"
DATE="$(date +"%y%m%d")"
TIME="$(date +"%H%M")"

# Get IDs of failed jobs
dx find jobs \
    --project project-GKK5xq0J7yj8yZZ863Jgg51x \
    --state failed \
    --brief \
    --name "${JOB_NAME}" \
    --tag "${TAG_01}" \
    -n 1000 >\
    "${JOB_ID_FILE}"

# Rerun the jobs
while read -r JOB_ID
    do
        dx run \
            --clone "${JOB_ID}" \
            --instance-type="${instance_type}" \
            --instance-count="${instances}" \
            --cost-limit="${cost_limit}" \
            --tag "${DATE}" \
            --tag "${TIME}" \
            --tag "${RERUN}" \
            -y \
            --brief
    done < "${JOB_ID_FILE}"

# Clean up
rm "${JOB_ID_FILE}"