.ONESHELL:
.PHONY: all test

SHELL=/bin/bash

vep_vcf = data/interim/cds_all_possible_snvs_vep.vcf.gz \
          data/interim/cds_all_possible_snvs_vep.vcf.gz.tbi \
		  data/manual/header_lines_vep_consequence.txt \

gnomad_tsv = data/manual/header_lines_gnomad.txt \
             data/interim/gnomad_v4.1_pass_snvs.tsv.gz \

all : $(vep_vcf) \
      $(gnomad_tsv) \

$(vep_vcf) &: data/interim/cds_all_possible_snvs_vep.tsv.gz \
              src/snv_annotation/vep_tsv_to_vcf.sh
	source activate vep
	bash src/snv_annotation/vep_tsv_to_vcf.sh
	# Takes about 5 minutes to run

$(gnomad_tsv) &: data/interim/gnomad_v4.1_pass_snvs.vcf.gz \
                 src/snv_annotation/gnomad_pass_snvs_to_tsv.sh
	source activate vep
	bash src/snv_annotation/gnomad_pass_snvs_to_tsv.sh
	
# # Merge annotations for all SNVs
# data/final/all_variants_merged_annotations.tsv : data/interim/nmd_annotations.tsv \
#                                                  data/interim/cds_all_possible_snvs.vcf \
# 												 data/interim/cds_trinucleotide_contexts.tsv \
# 												 data/interim/cds_all_possible_snvs_vep_tidy.tsv \
# 												 data/interim/gnomad_v4_pass_snvs.tsv \
# 												 data/interim/mutation_rate_by_context_methyl_tidy.tsv \
# 												 src/data/observed_variants.py
# 	python3 -m src.data.observed_variants
