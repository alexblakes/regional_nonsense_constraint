.ONESHELL:
.PHONY: all test

SHELL=/bin/bash

all : data/interim/gencode_v39_canonical_cds_seq.tsv \
      data/interim/cds_all_possible_snvs.vcf \
	  data/interim/cds_all_possible_snvs.vcf.gz \
	  data/interim/cds_trinucleotide_contexts.tsv \
	  data/interim/gnomad_v4.1_pass_snvs.vcf.gz \
	  data/interim/cds_all_possible_snvs_vep.vcf.gz \
	  data/logs/all_possible_snvs.log \

# Get FASTA sequences for CDS regions
data/interim/gencode_v39_canonical_cds_seq.tsv : data/interim/gencode_v39_canonical_cds.bed \
												 src/snvs/get_fasta.sh
	source activate vep
	bash src/snvs/get_fasta.sh

# Get all possible SNVs and trinucleotide contexts
data/interim/cds_all_possible_snvs.vcf \
data/interim/cds_trinucleotide_contexts.tsv &: data/interim/gencode_v39_canonical_cds_seq.tsv \
                                               src/snvs/coding_snvs.py
	source activate ukb
	python3 -m src.snvs.coding_snvs

# Get all possible SNVs (VCF format)
data/interim/cds_all_possible_snvs.vcf.gz : data/interim/gencode_v39_canonical_cds.bed \
                                            src/snvs/all_possible_snvs.sh
	source activate vep
	bash src/snvs/all_possible_snvs.sh
	# Takes about 30 minutes to run

# Extract SNVs from gnomAD
data/interim/gnomad_v4.1_pass_snvs.vcf.gz : data/interim/gencode_v39_canonical_cds.bed \
                                            src/snvs/extract_pass_variants.sh
	source activate vep
	bash src/snvs/extract_pass_variants.sh
	# Takes about 60 minutes to run

# Annotate all possible CDS SNVs with VEP
data/interim/cds_all_possible_snvs_vep.vcf.gz : data/interim/cds_all_possible_snvs.vcf.gz \
												src/snvs/vep_all_snvs_submit.sh \
												src/snvs/vep_all_snvs.sh \
									            src/snvs/vep_all_snvs_combine.sh
	source activate vep
	bash src/snvs/vep_all_snvs_submit.sh
	# Run from the CSF login node
	# This takes about 3 hours, if all jobs run immediately.

# Log some key statistics from the SNV annotation pipeline.
data/logs/all_possible_snvs.log : src/snvs/logging.sh
	bash src/snvs/logging.sh &

# # Tidy mutability data
# data/interim/mutation_rate_by_context_methyl_tidy.tsv : src/data/mutability_data.py
# 	python3 -m src.data.mutability_data

# # Merge annotations for all SNVs
# data/final/all_variants_merged_annotations.tsv : data/interim/nmd_annotations.tsv \
#                                                  data/interim/cds_all_possible_snvs.vcf \
# 												 data/interim/cds_trinucleotide_contexts.tsv \
# 												 data/interim/cds_all_possible_snvs_vep_tidy.tsv \
# 												 data/interim/gnomad_v4_pass_snvs.tsv \
# 												 data/interim/mutation_rate_by_context_methyl_tidy.tsv \
# 												 src/data/observed_variants.py
# 	python3 -m src.data.observed_variants
