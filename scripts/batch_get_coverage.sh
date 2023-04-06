#!/usr/bin/env bash

# This script submit batch jobs to extract coverage statistics from the UKB exomes VCFs

# Identifier variables for job tags
DATE="$(date +"%y%m%d")"
TIME="$(date +"%H%M")"

# Create local directories
rm -rf split_paths
mkdir split_paths

# Create UKB RAP directories
dx rm -rf outputs/gnomad_coverage/
dx mkdir -p outputs/gnomad_coverage/split_coverage 
dx mkdir outputs/gnomad_coverage/output_notebooks

# Get file paths
f=./gnomad_exome_file_paths.txt
if test -f "${f}"; then
    echo "${f} already exists."
else
    find "/mnt/project/Bulk/Exome sequences_Alternative exome processing/Exome variant call files (gnomAD) (VCFs)/" |\
    grep ".vcf.gz" |\
    grep -v ".tbi" |\
    sort >\
    gnomad_exome_file_paths.txt
fi

# Header file for testing purposes
head -n 50 gnomad_exome_file_paths.txt > test_file_paths.txt

# Create variable for the filename containing all VCF file paths
ALL_PATHS="test_file_paths.txt"

# Upload all VCF file paths to UKB RAP
dx upload --destination /outputs/gnomad_coverage/ "${ALL_PATHS}"

# Split ALL_PATHS into N VCFs per file
N=10
split -l "${N}" "${ALL_PATHS}" split_paths/split_paths_

# Upload split files to UKB RAP
dx upload -r --destination /outputs/gnomad_coverage/ split_paths

# Submit batch jobs
for LOCAL_FILE_PATH in split_paths/*;
do
    FILE_NAME="$(basename "${LOCAL_FILE_PATH}")"
    dx run dxjupyterlab_spark_cluster \
    -iin="${project}:/scripts/get_coverage.ipynb" \
    -iin="${project}:/outputs/gnomad_coverage/split_paths/${FILE_NAME}" \
    -icmd="papermill get_coverage.ipynb ${FILE_NAME}_out.ipynb" \
    -ifeature=HAIL-0.2.78 \
    --name="get_coverage_test_230406_${FILE_NAME}" \
    --instance-type=mem3_ssd1_v2_x2 \
    --instance-count=1 \
    --destination="${project}:/outputs/gnomad_coverage/output_notebooks/" \
    --priority="low" \
    --cost-limit 1.00 \
    --tag "${DATE}" \
    --tag "${TIME}" \
    --tag "test_01" \
    --tag "mem3_ssd1_v2_x2" \
    --tag "n_jobs=5" \
    --tag "n_vcfs=50" \
    -y
done

# Clean up
#rm -rf split_paths/ gnomad_exome_file_paths.txt test_file_paths.txt
