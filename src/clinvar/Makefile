.ONESHELL:
.PHONY: all

SHELL = bash

all : data/interim/clinvar_variants_selected.tsv \
      data/interim/clinvar_variants_selected.vcf \
	  data/interim/clinvar_variants_annotated.vcf \
	  data/interim/clinvar_variants_vep.vcf \
	  data/interim/clinvar_variants_vep_tidy.tsv \
	  data/interim/clinvar_variants_constraint.tsv \
	  data/plots/clinvar_acmg.png \

# Parse ClinVar summary text file
data/interim/clinvar_variants_selected.tsv : data/raw/variant_summary.txt \
                                             src/clinvar/filter_clinvar_variants.py
	source activate ukb
	python3 -m src.clinvar.filter_clinvar_variants

# Convert variants to VCF
data/interim/clinvar_variants_selected.vcf \
data/interim/clinvar_variants_annotated.vcf &: src/clinvar/to_vcf.sh \
											   src/clinvar/annotate_vcf.sh \
                                               data/interim/clinvar_variants_selected.tsv \
											   data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
	source activate bio
	bash src/clinvar/to_vcf.sh
	bash src/clinvar/annotate_vcf.sh
	source activate ukb

# Annotate with VEP
data/interim/clinvar_variants_vep.vcf : data/interim/clinvar_variants_annotated.vcf \
                                        src/clinvar/vep.sh
	source activate vep
	bash src/clinvar/vep.sh
	source activate ukb

# Tidy VEP output
data/interim/clinvar_variants_vep_tidy.tsv : data/interim/clinvar_variants_vep.vcf \
                                             data/final/nmd_annotations_simple.tsv.gz \
											 data/manual/header_lines_nmd_annotation.txt \
											 src/clinvar/tidy_vep.sh
	source activate bio
	bash src/clinvar/tidy_vep.sh
	source activate ukb

# Annotate with constraint
data/interim/clinvar_variants_constraint.tsv : data/interim/clinvar_variants_vep_tidy.tsv \
                                               data/final/regional_nonsense_constraint.tsv \
											   src/clinvar/annotate_with_constraint.py
	python3 -m src.clinvar.annotate_with_constraint

# Plot ACMG proportions
data/plots/clinvar_acmg.png : data/interim/clinvar_variants_constraint.tsv \
                              src/clinvar/plot_acmg.py
	python3 -m src.clinvar.plot_acmg