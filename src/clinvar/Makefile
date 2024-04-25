.ONESHELL:
.PHONY: all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

all : data/interim/clinvar_variants_selected.tsv \
      data/interim/clinvar_variants_selected.vcf \
	  data/interim/clinvar_variants_annotated.vcf \
	  data/final/nmd_annotations_simple.tsv.gz \
	  data/final/nmd_annotations_simple.tsv.gz.tbi \

# Parse ClinVar summary text file
data/interim/clinvar_variants_selected.tsv : data/raw/variant_summary.txt \
                                             src/clinvar/filter_clinvar_variants.py
	$(CONDA_ACTIVATE) ukb
	python3 -m src.data.filter_clinvar_variants

# Convert variants to VCF
data/interim/clinvar_variants_selected.vcf \
data/interim/clinvar_variants_annotated.vcf &: src/clinvar/to_vcf.sh \
                                               data/interim/clinvar_variants_selected.tsv \
											   data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
											   data/manual/clinvar_vcf_header_lines.txt
	$(CONDA_ACTIVATE) bio
	bash src/clinvar/to_vcf.sh
	bash src/clinvar/annotate_vcf.sh
	$(CONDA_ACTIVATE) ukb

# Create indexed NMD annotations
data/final/nmd_annotations_simple.tsv.gz \
data/final/nmd_annotations_simple.tsv.gz.tbi : src/clinvar/nmd_annotation.sh \
                                               data/interim/nmd_annotations.tsv
	$(CONDA_ACTIVATE) bio
	bash src/clinvar/nmd_annotation.sh
	$(CONDA_ACTIVATE) ukb


# Annotate with VEP
# $(CONDA_ACTIVATE) vep
# $(CONDA_ACTIVATE) ukb


# # VEP-annotate ClinVar variants
# data/interim/clinvar_variants_vep.tsv : data/interim/clinvar_variants_selected.vcf \
#                                         src/data/clinvar_vep.sh
# 	$(CONDA_ACTIVATE) vep
# 	bash src/data/clinvar_vep.sh
# 	$(CONDA_ACTIVATE) ukb

# # Tidy VEP-annotated ClinVar variants
# data/interim/clinvar_variants_vep_tidy.tsv : data/interim/clinvar_variants_vep.tsv \
#                                              src/data/clinvar_vep_tidy.py
# 	python3 -m src.data.clinvar_vep_tidy

# # Annotate ClinVar variants with regional constraint
# data/interim/clinvar_variants_lof_with_nmd_annotation.tsv : data/interim/clinvar_variants_vep_tidy.tsv \
#                                                             src/functional_clinical/clinvar_variants_in_constrained_regions.py
# 	python3 -m src.functional_clinical.clinvar_variants_in_constrained_regions
