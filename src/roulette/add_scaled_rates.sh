#!/usr/bin/env bash
source activate roulette
set -euo pipefail

# Run the Python script to scale mutation rates

VCF_DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs"
INPUT="data/interim/gnomad_v4.1_pass_snvs_roulette.tsv"
OUTPUT_DIR="data/interim/roulette"

python src/roulette/add_scaled_rates.py --vcf_dir $VCF_DIR --input $INPUT --output_dir $OUTPUT_DIR