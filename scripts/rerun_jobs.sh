#!/usr/bin/env bash

# This script reruns batch jobs which have failed

# Change these variables for each batch of jobs
JOB_NAME="*fail*original*230226*"
RERUN="rerun1"

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
    -n 1000 >\
    "${JOB_ID_FILE}"

# Rerun the jobs
while read -r JOB_ID
    do
        dx run \
            --clone "${JOB_ID}" \
            --tag "${DATE}" \
            --tag "${TIME}" \
            --tag "${RERUN}" \
            -y \
            --brief
    done < "${JOB_ID_FILE}"

# Clean up
rm "${JOB_ID_FILE}"