.ONESHELL:
.PHONY: all test

SHELL=/bin/bash

all : 

# Refactor gnomAD observed SNVs for Roulette input
data/interim/gnomad_v4.1_pass_snvs_roulette.tsv : data/interim/gnomad_v4.1_pass_snvs.vcf.gz \
                                                  src/roulette/tidy_gnomad_snvs.sh
	source activate vep
	bash src/roulette/tidy_gnomad_snvs.sh

# Scale Roulette mutation rates
test : src/roulette/add_scaled_rates.py
	source activate roulette
	python src/roulette/add_scaled_rates.py --vcf /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs --input data/interim/gnomad_v4.1_pass_snvs_roulette.tsv --output_dir data/interim/roulette