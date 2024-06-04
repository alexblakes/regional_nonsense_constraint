.ONESHELL:
.PHONY: all test

SHELL=/bin/bash

scaled = data/interim/roulette/1_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/2_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/3_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/4_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/5_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/6_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/7_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/8_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/9_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/10_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/11_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/12_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/13_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/14_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/15_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/16_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/17_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/18_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/19_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/20_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/21_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \
         data/interim/roulette/22_scaled_rate_v5.2_TFBS_correction_all.tsv.gz \

scaled_bgz = data/interim/roulette/1_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/2_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/3_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/4_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/5_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/6_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/7_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/8_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/9_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/10_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/11_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/12_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/13_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/14_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/15_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/16_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/17_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/18_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/19_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/20_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/21_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \
             data/interim/roulette/22_scaled_rate_v5.2_TFBS_correction_all.tsv.bgz \

raw_indices = /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/1_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/2_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/3_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/4_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/5_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/6_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/7_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/8_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/9_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/10_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/11_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/12_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/13_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/14_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/15_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/16_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/17_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/18_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/19_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/20_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/21_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \
              /mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs/22_rate_v5.2_TFBS_correction_all.vcf.bgz.tbi \

all : data/interim/gnomad_v4.1_pass_snvs_roulette.tsv \
      $(scaled) \
      $(scaled_bgz) \
      data/interim/roulette/scaled_cds_combined.tsv.gz \
      data/interim/roulette/scaled_cds_combined_sorted.tsv.gz \
      $(raw_indices) \
      data/interim/roulette/raw_cds_combined.vcf.gz \
      data/interim/roulette/raw_cds_combined.tsv.gz \
      data/logs/roulette_scores.log \

# Refactor gnomAD observed SNVs for Roulette input
data/interim/gnomad_v4.1_pass_snvs_roulette.tsv : data/interim/gnomad_v4.1_pass_snvs.vcf.gz \
                                                  src/roulette/tidy_gnomad_snvs.sh
	source activate vep
	bash src/roulette/tidy_gnomad_snvs.sh

# Scale Roulette mutation rates
$(scaled) &: data/interim/gnomad_v4.1_pass_snvs_roulette.tsv \
             src/roulette/submit_scaling.sh \
			 src/roulette/add_scaled_rates.sh \
             src/roulette/add_scaled_rates.py
	source activate roulette
	bash src/roulette/submit_scaling.sh
      # Run from the CSF
      # This job takes >24 hours to run.

# Tidy scaled mutation rates
$(scaled_bgz) &: src/roulette/scaled_drop_header.sh
	source activate vep
	bash src/roulette/scaled_drop_header.sh
	# Takes about 10 minutes to run

# Get the scaled Roulette scores in CDS sites
data/interim/roulette/scaled_cds_combined.tsv.gz : $(scaled_bgz) \
                                                   data/interim/gencode_v39_canonical_cds.bed \
                                                   src/roulette/scaled_concat.sh
	source activate vep
	bash src/roulette/scaled_concat.sh
	# Takes about 5 minutes to run

# Sort and index the scaled Roulette scores
data/interim/roulette/scaled_cds_combined_sorted.tsv.gz : data/interim/roulette/scaled_cds_combined.tsv.gz \
                                                          src/roulette/scaled_to_tsv.sh
	source activate vep
	bash src/roulette/scaled_to_tsv.sh
	# Takes about 20 minutes to run

# Tabix index the raw Roulette scores
$(raw_indices) &: src/roulette/raw_index_bgz_files.sh
	source activate vep
	bash src/roulette/raw_index_bgz_files.sh
	# Takes about 10 minutes to run

# Get raw Roulette scores in CDS sites
data/interim/roulette/raw_cds_combined.vcf.gz : data/interim/gencode_v39_canonical_cds_no_chr.bed \
                                                src/roulette/raw_concat.sh
	source activate vep
	bash src/roulette/raw_concat.sh
	# Takes about 10 minutes to run

# Convert raw Roulette scores to TSV
data/interim/roulette/raw_cds_combined.tsv.gz : data/interim/roulette/raw_cds_combined.vcf.gz \
                                                src/roulette/raw_to_tsv.sh
	source activate vep
	bash src/roulette/raw_to_tsv.sh
	# Takes about 10 minutes to run

# Log summary data about the processed Roulette data
data/logs/roulette_scores.log : data/interim/roulette/raw_cds_combined.tsv.gz \
                                data/interim/roulette/scaled_cds_combined_sorted.tsv.gz \
                                src/roulette/log.sh
	source activate vep
	bash src/roulette/log.sh
	# Takes about 5 minutes to run