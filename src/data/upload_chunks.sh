#!/usr/bin/env/ bash

# This script uploads the context_prepared.ht data in batches to the UKB RAP.
# Takes arguments:
# - $1 Directory to upload
# - $2 Number of files per chunk

# Get dx upload agent
wget \
  https://dnanexus-sdk.s3.amazonaws.com/dnanexus-upload-agent-1.5.33-linux.tar.gz -O - |\
  tar -xzf -

echo "DNAnexus upload agent ready."

# Get file paths for content of the hail table
find $1 -type f | sort > targets.txt

echo "List of targets ready."

# Create chunks to upload
rm -rf split/
mkdir split
split -l $2 targets.txt split/

echo "Split files ready."

# Iteratively upload files by chunk
for SPLIT_FILE in split/*
do
    while read TARGET;
    do
        DIR="$(dirname "${TARGET}")"
        dnanexus-upload-agent-1.5.33-linux/ua --do-not-compress -f "/data/${DIR}/" "${TARGET}" &
    done < "${SPLIT_FILE}"
    wait
    echo "Split file $(basename "${SPLIT_FILE}") done!"
done;

touch _UPLOAD_SUCCESS
dx upload --destination data/ _UPLOAD_SUCCESS