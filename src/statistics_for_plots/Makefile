.ONESHELL:
.PHONY: all

SHELL = bash

all : data/statistics/regions_cds_proportions.tsv \
	  data/statistics/constraint_upset_data.tsv \
      data/statistics/clinvar_vus_by_region.tsv \
	  data/statistics/clinvar_ptv_ascertainment.tsv \
	  data/statistics/clinvar_acmg_by_region_and_constraint.tsv \

# Get relative footprint of NMD regions
data/statistics/regions_cds_proportions.tsv : data/final/nmd_annotations_simple.tsv.gz \
                                              src/statistics_for_plots/regions_cds_proportions.py
	source activate ukb
	python3 -m src.statistics_for_plots.regions_cds_proportions

# Reformat constraint data for upset plot
data/statistics/constraint_upset_data.tsv : data/final/regional_nonsense_constraint.tsv \
                                            src/statistics_for_plots/constraint_upset.py
	source activate ukb
	python3 -m src.statistics_for_plots.constraint_upset

# Get proportion VUS in ClinVar by region
data/statistics/clinvar_vus_by_region.tsv : data/interim/clinvar_variants_vep_tidy.tsv \
                                            src/statistics_for_plots/clinvar_proportion_vus_by_region.py
	source activate ukb
	python3 -m src.statistics_for_plots.clinvar_proportion_vus_by_region

# Get PTV ascertainment in ClinVar by region
data/statistics/clinvar_ptv_ascertainment.tsv : data/interim/clinvar_variants_vep_tidy.tsv \
                                                data/statistics/regions_cds_proportions.tsv \
												src/statistics_for_plots/clinvar_ptv_ascertainment.py
	source activate ukb
	python3 -m src.statistics_for_plots.clinvar_ptv_ascertainment

# Get ACMG classification in ClinVar by region and constraint
data/statistics/clinvar_acmg_by_region_and_constraint.tsv : data/interim/clinvar_variants_constraint.tsv \
                                                            src/statistics_for_plots/clinvar_acmg_by_region_and_constraint.py
	source activate ukb
	python3 -m src.statistics_for_plots.clinvar_acmg_by_region_and_constraint
