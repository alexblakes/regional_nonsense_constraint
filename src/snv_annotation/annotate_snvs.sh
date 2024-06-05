#!/usr/bin/env bash
set -euo pipefail

# Annotate SNVs

VEP_VCF="data/interim/cds_all_possible_snvs_vep.vcf.gz"
VEP_TSV="data/interim/cds_all_possible_snvs_vep.tsv.gz"
NMD_TSV="data/final/nmd_annotations_simple.tsv.gz"
GNOMAD_TSV="data/interim/gnomad_v4.1_pass_snvs.tsv.gz"
ROULETTE_RAW_TSV="data/interim/roulette/raw_cds_combined.tsv.gz"
ROULETTE_SCALED_TSV="data/interim/roulette/scaled_cds_combined_sorted.tsv.gz"
HEADER_VEP="data/manual/header_lines_vep_consequence.txt"
HEADER_NMD="data/manual/header_lines_nmd_annotation.txt"
HEADER_GNOMAD="data/manual/header_lines_gnomad.txt"
HEADER_ROULETTE_RAW="data/manual/header_lines_roulette_raw_bgz.txt"
HEADER_ROULETTE_SCALED="data/manual/header_lines_roulette_scaled.txt"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# The -a file must be in tsv format in order for bcftools annotate to work in a stream.

# Annotate SNVs
bcftools annotate -a $VEP_TSV -c CHROM,POS,REF,ALT,~ID,Consequence -h $HEADER_VEP -Ou $VEP_VCF \
| bcftools annotate -a $NMD_TSV -c CHROM,POS,~ID,REGION -h $HEADER_NMD -Ou \
| bcftools annotate -a $GNOMAD_TSV -c CHROM,POS,REF,ALT,AC,AN,AF -h $HEADER_GNOMAD \
| grep -v "^##" | less -S

# # Merge annotations for all SNVs
# data/final/all_variants_merged_annotations.tsv : ***data/interim/nmd_annotations.tsv \
#                                                  ***data/interim/cds_all_possible_snvs.vcf \
# 												 ---data/interim/cds_trinucleotide_contexts.tsv \
# 												 ***data/interim/cds_all_possible_snvs_vep_tidy.tsv \
# 												 ***data/interim/gnomad_v4_pass_snvs.tsv \
# 												 data/interim/mutation_rate_by_context_methyl_tidy.tsv \
# 												 ---src/data/observed_variants.py
# 	python3 -m src.data.observed_variants
