.ONESHELL:
.PHONY: all

SHELL = /bin/bash

all : data/interim/gencode_v39_canonical_cds_seq.tsv \
      data/interim/cds_all_possible_snvs.vcf \

# Get FASTA sequences for CDS regions
data/interim/gencode_v39_canonical_cds_seq.tsv : data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
                                                 data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai \
												 data/interim/gencode_v39_canonical_cds.bed \
												 src/snvs/get_fasta.sh
	source activate bio
	bash src/snvs/get_fasta.sh

# Get all possible SNVs and trinucleotide contexts
data/interim/cds_all_possible_snvs.vcf :      data/interim/gencode_v39_canonical_cds_seq.tsv \
                                              src/snvs/coding_snvs.py
	source activate ukb
	python3 -m src.snvs.coding_snvs

# data/interim/cds_trinucleotide_contexts.tsv : data/interim/gencode_v39_canonical_cds_seq.tsv \
#                                               src/data/coding_snvs.py
# 	python3 -m src.data.coding_snvs

# # Extract SNVs from gnomAD
# data/interim/gnomad_v4_pass_snvs.tsv : data/interim/gencode_v39_canonical_cds.bed \
#                                        src/data/batch_extract_pass_variants.sh \
# 									   src/data/extract_pass_variants.sh \
# 									   src/data/combine_pass_snvs.sh
# 	bash src/data/batch_extract_pass_variants.sh
# 	qsub src/data/combine_pass_snvs.sh

# # Annotate all possible CDS SNVs with VEP
# # This script runs over ~ 10 hours
# data/interim/vep_all_snvs/out_29.tsv : data/interim/cds_all_possible_snvs.vcf \
# 									   src/data/batch_vep_all_snvs.sh \
# 									   src/data/vep_all_snvs.sh
# 	bash src/data/batch_vep_all_snvs.sh
# 	qsub -hold_jid "vep_snvs_*"

# # Combine and tidy the split VEP outputs
# data/interim/cds_all_possible_snvs_vep_tidy.tsv : data/interim/vep_all_snvs/out_29.tsv \
#                                                   data/interim/transcript_ids.tsv \
#                                                   src/data/vep_all_snvs_tidy.sh \
# 												  src/data/vep_all_snvs_tidy_log.py
# 	bash src/data/vep_all_snvs_tidy.sh
# 	python3 -m src.data.vep_all_snvs_tidy_log

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
