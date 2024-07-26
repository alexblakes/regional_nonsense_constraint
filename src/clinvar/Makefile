.ONESHELL:
.PHONY: all

SHELL = bash

vus_plots = data/plots/clinvar/clinvar_vus_by_region.png \
            data/plots/clinvar/clinvar_vus_by_region.svg \

ptv_ascertainment_plots = data/plots/clinvar/clinvar_ptvs_ascertainment.png \
                          data/plots/clinvar/clinvar_ptvs_ascertainment.svg \

acmg_plots = data/plots/clinvar/clinvar_acmg.png \
             data/plots/clinvar/clinvar_acmg.svg \

all : data/interim/clinvar_variants_selected.tsv \
      data/interim/clinvar_variants_selected.vcf \
	  data/interim/clinvar_variants_annotated.vcf \
	  data/interim/clinvar_variants_vep.vcf \
	  data/interim/clinvar_variants_vep_tidy.tsv \
	  data/interim/clinvar_variants_constraint.tsv \
	  data/statistics/clinvar_vus_by_region.tsv \
	  data/statistics/clinvar_vus_by_region.tsv \
	  $(vus_plots) \
	  data/statistics/clinvar_ptv_ascertainment.tsv \
	  $(ptv_ascertainment_plots) \
	  $(acmg_plots) \

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
	source activate vep
	bash src/clinvar/to_vcf.sh
	bash src/clinvar/annotate_vcf.sh

# Annotate with VEP
data/interim/clinvar_variants_vep.vcf : data/interim/clinvar_variants_annotated.vcf \
                                        src/clinvar/vep.sh
	source activate vep
	bash src/clinvar/vep.sh
	# Takes ~30 minutes to run

# Tidy VEP output
data/interim/clinvar_variants_vep_tidy.tsv : data/interim/clinvar_variants_vep.vcf \
                                             data/final/nmd_annotations_simple.tsv.gz \
											 data/manual/header_lines_nmd_annotation.txt \
											 src/clinvar/tidy_vep.sh
	source activate vep
	bash src/clinvar/tidy_vep.sh

# Annotate with constraint
data/interim/clinvar_variants_constraint.tsv : data/interim/clinvar_variants_vep_tidy.tsv \
                                               data/final/regional_nonsense_constraint.tsv \
											   src/clinvar/annotate_with_constraint.py
	source activate ukb
	python3 -m src.clinvar.annotate_with_constraint