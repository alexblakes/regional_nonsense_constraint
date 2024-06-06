.ONESHELL:
.PHONY: all force

SHELL=/bin/bash

all : data/interim/gencode_v39_canonical_cds_seq.tsv \
      data/interim/cds_all_possible_snvs.vcf \
	  data/interim/cds_all_possible_snvs.vcf.gz \
	  data/interim/cds_trinucleotide_contexts.tsv \
	  data/interim/gnomad_v4.1_pass_snvs.vcf.gz \
	  data/interim/cds_all_possible_snvs_vep.tsv.gz \
	  data/logs/all_possible_snvs.log \

force: 

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
data/interim/cds_all_possible_snvs_vep.tsv.gz : data/interim/cds_all_possible_snvs.vcf.gz \
												src/snvs/vep_all_snvs_submit.sh \
												src/snvs/vep_all_snvs.sh \
									            src/snvs/vep_all_snvs_combine.sh
	source activate vep
	bash src/snvs/vep_all_snvs_submit.sh
	# Run from the CSF login node
	# Higher --buffer-size or --forks risk out-of-memory errors
	# This takes about 5 hours, if all jobs run immediately

# Log some key statistics from the SNV annotation pipeline.
data/logs/all_possible_snvs.log : src/snvs/logging_submit.sh \
                                  src/snvs/logging.sh \
                                  data/interim/cds_all_possible_snvs_vep.tsv.gz \
								  data/interim/cds_all_possible_snvs.vcf.gz
	source activate vep
	bash src/snvs/logging_submit.sh