.ONESHELL:
.PHONY: all

SHELL = bash

dominant_split_regions = clinvar_variants_in_dominant_genes_distal_nmd.tsv.gz \
                         clinvar_variants_in_dominant_genes_long_exon.tsv.gz \
                         clinvar_variants_in_dominant_genes_nmd_target.tsv.gz \
                         clinvar_variants_in_dominant_genes_start_proximal.tsv.gz

recessive_split_regions = clinvar_variants_in_recessive_genes_distal_nmd.tsv.gz \
                         clinvar_variants_in_recessive_genes_long_exon.tsv.gz \
                         clinvar_variants_in_recessive_genes_nmd_target.tsv.gz \
                         clinvar_variants_in_recessive_genes_start_proximal.tsv.gz

all : data/interim/clinvar_parsed.tsv \
      data/interim/mane_select_transcript_ids.tsv \
	  data/interim/clinvar_variants_selected.tsv \
      data/interim/clinvar_variants_selected.vcf.gz \
	  data/interim/clinvar_variants_annotated.vcf.gz \
	  data/interim/clinvar_variants_vep.tsv.gz \
	  data/interim/clinvar_variants_constraint.tsv.gz \
	  data/interim/clinvar_variants_in_dominant_genes.tsv.gz \
	  data/interim/clinvar_variants_in_recessive_genes.tsv.gz \
	  data/interim/clinvar_unexpected_variants_by_acmg_and_constraint.tsv \
	  data/interim/clinvar_unexpected_variants_by_acmg_and_constraint_sample.tsv

# Parse ClinVar summary text file
data/interim/clinvar_parsed.tsv : data/raw/variant_summary.txt \
                                  src/clinvar/filter_clinvar_variants.py
	source activate ukb
	python3 -m src.clinvar.filter_clinvar_variants


# Get MANE Select RefSeq and Ensembl transcript IDs
data/interim/mane_select_transcript_ids.tsv : data/raw/MANE.GRCh38.v0.95.summary.txt.gz \
                                              src/clinvar/get_mane_transcript_ids.py
	source activate ukb
	python3 -m src.clinvar.get_mane_transcript_ids


# Filter for ClinVar variants in MANE transcripts
data/interim/clinvar_variants_selected.tsv : data/interim/mane_select_transcript_ids.tsv \
                                             data/interim/clinvar_parsed.tsv \
											 src/clinvar/clinvar_variants_in_mane_transcripts.py
	source activate ukb
	python3 -m src.clinvar.clinvar_variants_in_mane_transcripts


# Convert variants to VCF
data/interim/clinvar_variants_selected.vcf.gz \
data/interim/clinvar_variants_annotated.vcf.gz &: src/clinvar/to_vcf.sh \
											      src/clinvar/annotate_vcf.sh \
                                                  data/interim/clinvar_variants_selected.tsv \
											      data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
	source activate vep
	bash src/clinvar/to_vcf.sh
	bash src/clinvar/annotate_vcf.sh


# Annotate with VEP
data/interim/clinvar_variants_vep.tsv.gz : data/interim/clinvar_variants_annotated.vcf.gz \
                                           data/final/nmd_annotations_simple.tsv.gz \
										   src/clinvar/vep.sh
	source activate vep
	bash src/clinvar/vep.sh
	# Takes ~30 minutes to run


# Annotate with constraint
data/interim/clinvar_variants_constraint.tsv.gz : data/interim/clinvar_variants_vep.tsv.gz \
                                               data/final/regional_nonsense_constraint.tsv \
											   src/clinvar/annotate_with_constraint.py
	source activate ukb
	python3 -m src.clinvar.annotate_with_constraint

# Filter for variants in morbid AD / AR genes
data/interim/clinvar_variants_in_dominant_genes.tsv.gz \
data/interim/clinvar_variants_in_recessive_genes.tsv.gz &: data/interim/clinvar_variants_constraint.tsv.gz \
                                                          data/interim/genemap2_simple.tsv \
														  src/clinvar/filter_for_morbid_genes.py
	source activate ukb
	python3 -m src.clinvar.filter_for_morbid_genes

# Split by region for dominant genes
$(dominant_split_regions) &: data/interim/clinvar_variants_in_dominant_genes.tsv.gz \
                             src/clinvar/split_by_region.py
	source activate ukb
	python3 -m src.clinvar.split_by_region data/interim/clinvar_variants_in_dominant_genes.tsv.gz

# Split by region for recessive genes
$(recessive_split_regions) &: data/interim/clinvar_variants_in_recessive_genes.tsv.gz \
                             src/clinvar/split_by_region.py
	source activate ukb
	python3 -m src.clinvar.split_by_region data/interim/clinvar_variants_in_recessive_genes.tsv.gz

# Get unexpected variants list
data/interim/clinvar_unexpected_variants_by_acmg_and_constraint.tsv \
data/interim/clinvar_unexpected_variants_by_acmg_and_constraint_sample.tsv &: data/interim/clinvar_variants_in_dominant_genes.tsv.gz \
                                                                              src/clinvar/list_unexpected_variants.py
	source activate ukb
	python3 -m src.clinvar.list_unexpected_variants

