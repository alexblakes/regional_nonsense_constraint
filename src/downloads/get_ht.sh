#!/usr/bin/env bash

# Script to get gnomAD non-coding constraint hail table.

# Install gsutil
conda install -c conda-forge gsutil -y

# Download hail table to worker
gsutil -m cp -r \
  "gs://gnomad-nc-constraint-v31-paper/context_prepared.ht" \
  .

