.ONESHELL:
.PHONY: all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate


# Parse ClinVar summary text file
data/interim/clinvar_variants_selected.tsv \
data/interim/clinvar_variants_selected.vcf : data/raw/variant_summary.txt \
                                             src/data/clinvar_variants.py
	python3 -m src.data.clinvar_variants

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
