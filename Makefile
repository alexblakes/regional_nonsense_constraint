.ONESHELL:
.PHONY: downloads \
        regions \
		snvs \
		roulette \
		snv_annotation \
		constraint \
		maps \
		cadd \
		clinvar \
		orthogonal_metrics \
		notebooks \
		all \
		# go \
		# gene_enrichment \

SHELL = /bin/bash

downloads: 
	make -f data/raw/Makefile all

regions :
	make -f src/regions/Makefile all

snvs :
	make -f src/snvs/Makefile all

roulette :
	make -f src/roulette/Makefile all

snv_annotation :
	make -f src/snv_annotation/Makefile all

constraint :
	make -f src/constraint/Makefile all

maps :
	make -f src/maps/Makefile all

cadd : 
	make -f src/cadd/Makefile all

clinvar :
	make -f src/clinvar/Makefile all

orthogonal_metrics :
	make -f src/orthogonal_metrics/Makefile all


statistics : 
	papermill src/statistics_for_plots/clinvar_ascertainment.ipynb src/statistics_for_plots/clinvar_ascertainment.ipynb
	papermill src/statistics_for_plots/maps.ipynb src/statistics_for_plots/maps.ipynb
	papermill src/statistics_for_plots/oe.ipynb src/statistics_for_plots/oe.ipynb
	papermill src/statistics_for_plots/z_loeuf.ipynb src/statistics_for_plots/z_loeuf.ipynb
	papermill src/statistics_for_plots/z_distributions.ipynb src/statistics_for_plots/z_distributions.ipynb
	papermill src/statistics_for_plots/upset.ipynb src/statistics_for_plots/upset.ipynb

figures :
	papermill notebooks/figures/fig_01.ipynb notebooks/figures/fig_01.ipynb

notebooks :
	# Expectation model choices
	for N in 0 10 20 30 ; do \
		papermill notebooks/01_expectation_model_choices.ipynb \
		notebooks/01_expectation_model_choices.ipynb \
		-p coverage $$N ; \
	done

# gene_enrichment :
# 	# Run this from the CSF; it requires API access to gProfiler
# 	make -f src/gene_enrichment/Makefile all

all : downloads \
	  regions \
	  snvs \
	  roulette \
	  snv_annotation \
	  constraint \
	  maps \
	  cadd \
	  clinvar \
	  orthogonal_metrics \
	  statistics \
	  figures \
      # gene_enrichment \

### ORTHOGONAL METRICS


### OMIM

# Parse genemap2.txt data
data/interim/genemap2_parsed.tsv : data/raw/genemap2.txt \
                                   src/functional_clinical/omim_parse_genemap.py
	python3 -m src.functional_clinical.omim_parse_genemap

# Simplify genemap2.txt data
data/interim/genemap2_simple.tsv : data/interim/genemap2_parsed.tsv \
                                   src/functional_clinical/omim_simplify_genemap.py
	python3 -m src.functional_clinical.omim_simplify_genemap

### CLINVAR

# Parse ClinVar summary text file
data/interim/clinvar_variants_selected.tsv \
data/interim/clinvar_variants_selected.vcf : data/raw/variant_summary.txt \
                                             src/data/clinvar_variants.py
	python3 -m src.data.clinvar_variants

# VEP-annotate ClinVar variants
data/interim/clinvar_variants_vep.tsv : data/interim/clinvar_variants_selected.vcf \
                                        src/data/clinvar_vep.sh
	$(CONDA_ACTIVATE) vep
	bash src/data/clinvar_vep.sh
	$(CONDA_ACTIVATE) ukb

# Tidy VEP-annotated ClinVar variants
data/interim/clinvar_variants_vep_tidy.tsv : data/interim/clinvar_variants_vep.tsv \
                                             src/data/clinvar_vep_tidy.py
	python3 -m src.data.clinvar_vep_tidy

# Annotate ClinVar variants with regional constraint
data/interim/clinvar_variants_lof_with_nmd_annotation.tsv : data/interim/clinvar_variants_vep_tidy.tsv \
                                                            src/functional_clinical/clinvar_variants_in_constrained_regions.py
	python3 -m src.functional_clinical.clinvar_variants_in_constrained_regions

# Statistics


# Next 