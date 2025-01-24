.ONESHELL:
.PHONY: all

SHELL = bash

all : data/statistics/regions_cds_proportions.tsv \
      data/statistics/region_counts.tsv \
	  data/statistics/constraint_upset_data.tsv \
      data/statistics/clinvar_vus_by_region.tsv \
	  data/statistics/clinvar_ptv_ascertainment.tsv \
	  data/statistics/clinvar_acmg_by_region_and_constraint.tsv \
	  data/statistics/clinvar_by_csq_constraint.tsv \
	  data/statistics/clinvar_by_csq_constraint_morbid_dominant.tsv \
	  data/statistics/orthogonal_metrics.tsv.gz \
	  data/statistics/orthogonal_metrics_cadd.tsv.gz \
	  data/statistics/maps.tsv \

# Get relative footprint of NMD regions
data/statistics/regions_cds_proportions.tsv : data/final/nmd_annotations_simple.tsv.gz \
                                              src/statistics_for_plots/regions_cds_proportions.py
	source activate ukb
	python3 -m src.statistics_for_plots.regions_cds_proportions

# Count transcripts with a given NMD region
data/statistics/region_counts.tsv : data/final/nmd_annotations_simple.tsv.gz \
                                    src/statistics_for_plots/region_counts.py
	source activate ukb
	python3 -m src.statistics_for_plots.region_counts

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

# Get proportion constrained for ClinVar variants split by consequence and ACMG class.
data/statistics/clinvar_by_csq_constraint.tsv : data/interim/clinvar_variants_constraint.tsv.gz \
                                                src/statistics_for_plots/clinvar_by_csq_constraint.py
	source activate ukb
	python3 -m src.statistics_for_plots.clinvar_by_csq_constraint data/interim/clinvar_variants_constraint.tsv.gz data/statistics/clinvar_by_csq_constraint.tsv

# Get proportion constrained for ClinVar variants split by consequence and ACMG class, for variants in morbid dominant genes only
data/statistics/clinvar_by_csq_constraint_morbid_dominant.tsv : data/interim/clinvar_variants_constraint_in_morbid_dominant_genes.tsv.gz \
                                                                src/statistics_for_plots/clinvar_by_csq_constraint.py
	source activate ukb
	python3 -m src.statistics_for_plots.clinvar_by_csq_constraint data/interim/clinvar_variants_constraint_in_morbid_dominant_genes.tsv.gz data/statistics/clinvar_by_csq_constraint_morbid_dominant.tsv

# Get phyloP, AlphaMissense, and pext scores
data/statistics/orthogonal_metrics.tsv.gz : data/interim/cds_sites_phylop_pext_missense.tsv \
                                            src/statistics_for_plots/orthogonal_metrics.py
	source activate ukb
	python3 -m src.statistics_for_plots.orthogonal_metrics

# Get CADD scores
data/statistics/orthogonal_metrics_cadd.tsv.gz : data/interim/cadd_scores_coding_annotated.tsv \
                                                 src/statistics_for_plots/orthogonal_metrics_cadd.py
	source activate ukb
	python3 -m src.statistics_for_plots.orthogonal_metrics_cadd

# MAPS
data/statistics/maps.tsv : data/interim/maps.tsv \
                           src/statistics_for_plots/maps.py
	source activate ukb
	python3 -m src.statistics_for_plots.maps