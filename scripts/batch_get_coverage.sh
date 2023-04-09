#!/usr/bin/env bash

# This script submit batch jobs to extract coverage statistics from the UKB exomes VCFs

# Variables for the job submission script
run_type="production" # "test" or "production"
tag_name="original_01"
n_vcfs=80
n_jobs=1
files_in_test=1 # n_vcfs * n_jobs
instance_type="mem3_ssd1_v2_x2"
instances=8
priority="low"
DATE="$(date +"%y%m%d")"
TIME="$(date +"%H%M")"

# Do the work in a temporary directory
rm -rf tmp
mkdir tmp
cd tmp

# Create local directories
rm -rf split_paths
mkdir split_paths

# Create UKB RAP directories
dx rm -rf outputs/gnomad_coverage/
dx mkdir -p outputs/gnomad_coverage/split_coverage 
dx mkdir outputs/gnomad_coverage/output_notebooks

# Get file paths
FILE_PATHS="gnomad_exome_file_paths.txt"

if test -f "${FILE_PATHS}"; then
    echo "${FILE_PATHS} already exists."
else
    dx_path="Bulk/Exome sequences_Alternative exome processing/Exome variant call files (gnomAD) (VCFs)/"
    prefix="/mnt/project/"
    
    dx find data --path "${dx_path}" |\
    cut -d \/ -f 5 | cut -d ' ' -f 1 |\
    while read line; do echo "${prefix}${dx_path}${line}"; done |\
    grep ".vcf.gz" |\
    grep -v ".tbi" |\
    sort >\
    ${FILE_PATHS}
fi

# In testing, only run on a few VCFs:
if test "${run_type}" == "test"; then
    head -n "${files_in_test}" "${FILE_PATHS}" | sponge "${FILE_PATHS}"
fi

# Split FILE_PATHS into n_vcfs per file
split -l "${n_vcfs}" "${FILE_PATHS}" split_paths/split_paths_

# Upload VCF file paths and split files to UKB RAP
dx upload --destination /outputs/gnomad_coverage/ "${FILE_PATHS}" --brief
dx upload -r --destination /outputs/gnomad_coverage/ split_paths --brief

# Submit batch jobs
for LOCAL_FILE_PATH in split_paths/*;
do
    FILE_NAME="$(basename "${LOCAL_FILE_PATH}")"
    dx run dxjupyterlab_spark_cluster \
    -iin="${project}:/scripts/get_coverage.ipynb" \
    -iin="${project}:/outputs/gnomad_coverage/split_paths/${FILE_NAME}" \
    -icmd="papermill get_coverage.ipynb ${FILE_NAME}_out.ipynb" \
    -ifeature=HAIL-0.2.78 \
    --name="get_coverage_${run_type}_${DATE}_${FILE_NAME}" \
    --instance-type="${instance_type}" \
    --instance-count="${instances}" \
    --destination="${project}:/outputs/gnomad_coverage/output_notebooks/" \
    --priority="${priority}" \
    --cost-limit 0.40 \
    --tag "${DATE}" \
    --tag "${TIME}" \
    --tag "${tag_name}" \
    --tag "${instance_type}" \
    --tag "instances=${instances}" \
    --tag "n_jobs=${n_jobs}" \
    --tag "n_vcfs=${n_vcfs}" \
    -y \
    --brief
done

# Clean up
cd ..
rm -r tmp