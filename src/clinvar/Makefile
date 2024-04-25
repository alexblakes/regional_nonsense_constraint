.ONESHELL:
.PHONY: all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

all : data/interim/clinvar_variants_selected.tsv \
      data/interim/clinvar_variants_selected.vcf \
	  data/interim/clinvar_variants_annotated.vcf \
	  data/final/nmd_annotations_simple.tsv.gz \
	  data/final/nmd_annotations_simple.tsv.gz.tbi \
	  data/interim/clinvar_variants_vep.vcf \
	  data/interim/clinvar_variants_vep_tidy.tsv \
	  data/interim/clinvar_variants_constraint.tsv \

# Parse ClinVar summary text file
data/interim/clinvar_variants_selected.tsv : data/raw/variant_summary.txt \
                                             src/clinvar/filter_clinvar_variants.py
	$(CONDA_ACTIVATE) ukb
	python3 -m src.data.filter_clinvar_variants

# Convert variants to VCF
data/interim/clinvar_variants_selected.vcf \
data/interim/clinvar_variants_annotated.vcf &: src/clinvar/to_vcf.sh \
											   src/clinvar/annotate_vcf.sh \
                                               data/interim/clinvar_variants_selected.tsv \
											   data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
	$(CONDA_ACTIVATE) bio
	bash src/clinvar/to_vcf.sh
	bash src/clinvar/annotate_vcf.sh
	$(CONDA_ACTIVATE) ukb

# Annotate with VEP
data/interim/clinvar_variants_vep.vcf : data/interim/clinvar_variants_annotated.vcf \
                                        src/clinvar/vep.sh
	$(CONDA_ACTIVATE) vep
	bash src/clinvar/vep.sh
	$(CONDA_ACTIVATE) ukb

# Create indexed NMD annotations
data/final/nmd_annotations_simple.tsv.gz \
data/final/nmd_annotations_simple.tsv.gz.tbi : src/clinvar/nmd_annotation.sh \
                                               data/interim/nmd_annotations.tsv
	$(CONDA_ACTIVATE) bio
	bash src/clinvar/nmd_annotation.sh
	$(CONDA_ACTIVATE) ukb

# Tidy VEP output
data/interim/clinvar_variants_vep_tidy.tsv : data/interim/clinvar_variants_vep.vcf \
                                             data/final/nmd_annotations_simple.tsv.gz \
											 data/manual/header_lines_nmd_annotation.txt \
											 src/clinvar/tidy_vep.sh
	$(CONDA_ACTIVATE) bio
	bash src/clinvar/tidy_vep.sh
	$(CONDA_ACTIVATE) ukb

# Annotate with constraint
data/interim/clinvar_variants_constraint.tsv : data/interim/clinvar_variants_vep_tidy.tsv \
                                               data/final/regional_nonsense_constraint.tsv \
											   src/clinvar/annotate_with_constraint.py
	python3 -m src.clinvar.annotate_with_constraint
