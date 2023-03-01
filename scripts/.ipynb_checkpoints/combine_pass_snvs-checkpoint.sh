#!/usr/bin/env bash

# This script downloads and combines all the UKB variants passing gnomAD filters

# Get filepaths readable to dx download
find "/mnt/project/outputs/gnomad_pass_variants/split_pass_variants/" -type f |\
sed 's/^.\{13\}//' |\
sort >\
pass_snps_file_paths_dx.txt

# Use dx download to download the files locally
# head -n 100 pass_snps_file_paths_dx.txt > test_paths.txt
rm -rf pass_snv_files/
mkdir -p pass_snv_files/
# Takes ~ 1 hour to run
xargs -a pass_snps_file_paths_dx.txt -n 10 -P 0 -I % dx download -o pass_snv_files/ %


# Concatenate the individual files containing PASS SNVs
# The number of files is too great for cat, so passing to xargs is required
find pass_snv_files/ -type f -name *.tsv | xargs -n 1000 -I % cat % > all_pass_snvs.txt

# Upload the combined data
dx rm -f outputs/gnomad_pass_variants/all_pass_snvs.txt
dx upload --destination outputs/gnomad_pass_variants/ all_pass_snvs.txt