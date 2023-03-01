#!/usr/bin/env bash

# This script submits batch jobs to extract variants from UKB exomes

# Download the DNAnexus Upload Agent
wget \
  https://dnanexus-sdk.s3.amazonaws.com/dnanexus-upload-agent-1.5.33-linux.tar.gz -O - |\
  tar -xzf -

# Create local directories
rm -rf split_paths
mkdir split_paths

# Create directories in UKB RAP
dx rm -rf outputs/gnomad_pass_variants/
dx mkdir -p outputs/gnomad_pass_variants/split_pass_variants

# Get file paths
find "/mnt/project/Bulk/Exome sequences_Alternative exome processing/Exome variant call files (gnomAD) (VCFs)/" |\
grep ".vcf.gz" |\
sed 's/^.\{12\}/project-GKK5xq0J7yj8yZZ863Jgg51x\:/' |\
sort >\
gnomad_exome_file_paths.txt

## Header file for testing purposes
head -n 8 gnomad_exome_file_paths.txt > test_file_paths.txt

ALL_PATHS="gnomad_exome_file_paths.txt"

dx upload --destination /outputs/gnomad_pass_variants/ "${ALL_PATHS}"

# Split ALL_PATHS into N/2 VCFs per file (NB ALL_PATHS includes .tbi indexes)
N=200
split -l "${N}" "${ALL_PATHS}" split_paths/split_paths_

# Upload split files to UKB RAP using the Upload Agent
dnanexus-upload-agent-*-linux/ua --do-not-compress -f /outputs/gnomad_pass_variants/ split_paths

# Submit batch jobs
for LOCAL_FILE_PATH in split_paths/*;
do
    FILE_NAME="$(basename "${LOCAL_FILE_PATH}")"
    dx run swiss-army-knife \
    -iin="${project}:/scripts/extract_pass_variants.sh" \
    -iin="${project}:/outputs/gnomad_pass_variants/split_paths/${FILE_NAME}" \
    -iin="${project}:/data/gencode_v39_canonical_cds_chr.bed" \
    -icmd="bash extract_pass_variants.sh" \
    --name="extract_pass_variants_original_230225 ${FILE_NAME}" \
    --instance-type="mem3_ssd1_v2_x8" \
    --destination="${project}:/outputs/gnomad_pass_variants/split_pass_variants/" \
    --priority="low" \
    --cost-limit 0.40 \
    --tag "original" \
    --tag "parallel" \
    --tag "mem3_ssd1_v2_x8" \
    --tag "dx_download" \
    -y
    sleep 1
done