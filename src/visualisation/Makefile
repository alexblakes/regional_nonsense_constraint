.ONESHELL:
.PHONY: all

SHELL = bash

all : data/plots/regions/regions_cds_proportions.png \
	  data/plots/constraint/upset_plot.png \
      data/plots/clinvar/clinvar_vus_by_region.png \
	  data/plots/clinvar/clinvar_ptvs_ascertainment.png \
	  data/plots/clinvar/clinvar_acmg_by_region_and_constraint.png \
	  data/plots/orthogonal_metrics/phylop_alphamis_pext.png \
	  data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_dominant.png \
	  data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_recessive.png \
	  data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_dominant_distal.png \
	  data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_recessive_distal.png \
	  data/plots/orthogonal_metrics/cadd.png \
	  data/plots/maps/maps.png \

# Plot proportion of CDS per region
data/plots/regions/regions_cds_proportions.png : data/statistics/regions_cds_proportions.tsv \
                                                 src/visualisation/regions_cds_proportions.py
	source activate ukb
	python3 -m src.visualisation.regions_cds_proportions

# Constraint upset plot
data/plots/constraint/upset_plot.png : data/statistics/constraint_upset_data.tsv \
                                       src/visualisation/constraint_upset_plot.py 
	source activate ukb
	python3 -m src.visualisation.constraint_upset_plot

# Plot proportion VUS per region
data/plots/clinvar/clinvar_vus_by_region.png : data/statistics/clinvar_vus_by_region.tsv \
                                               src/visualisation/clinvar_proportion_vus_by_region.py
	source activate ukb
	python3 -m src.visualisation.clinvar_proportion_vus_by_region

# Plot proportion of ClinVar variants in constrained regions, by consequence and ACMG classification
data/plots/clinvar/clinvar_acmg_by_csq_and_constraint.png : data/statistics/clinvar_by_csq_constraint.tsv \
                                                            src/visualisation/clinvar_by_csq_constraint.py
	source activate ukb
	python3 -m src.visualisation.clinvar_by_csq_constraint data/statistics/clinvar_by_csq_constraint.tsv data/plots/clinvar/clinvar_acmg_by_csq_and_constraint.png

# Plot proportion of ClinVar variants in constrained regions, by consequence and ACMG classification, for dominant disease genes only
data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_dominant.png : data/statistics/clinvar_by_csq_constraint_morbid_dominant.tsv \
                                                                     src/visualisation/clinvar_by_csq_constraint.py
	source activate ukb
	python3 -m src.visualisation.clinvar_by_csq_constraint data/statistics/clinvar_by_csq_constraint_morbid_dominant.tsv data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_dominant.png

# Plot proportion of ClinVar variants in constrained regions, by consequence and ACMG classification, for dominant disease genes only
data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_recessive.png : data/statistics/clinvar_by_csq_constraint_morbid_recessive.tsv \
                                                                      src/visualisation/clinvar_by_csq_constraint.py
	source activate ukb
	python3 -m src.visualisation.clinvar_by_csq_constraint data/statistics/clinvar_by_csq_constraint_morbid_recessive.tsv data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_recessive.png

# Plot proportion of ClinVar variants in constrained regions, by consequence and ACMG classification, for distal regions in dominant disease genes only
data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_dominant_distal.png : data/statistics/clinvar_by_csq_constraint_morbid_dominant_distal.tsv \
                                                                            src/visualisation/clinvar_by_csq_constraint.py
	source activate ukb
	python3 -m src.visualisation.clinvar_by_csq_constraint data/statistics/clinvar_by_csq_constraint_morbid_dominant_distal.tsv data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_dominant_distal.png

# Plot proportion of ClinVar variants in constrained regions, by consequence and ACMG classification, for distal regions in recessive disease genes only
data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_recessive_distal.png : data/statistics/clinvar_by_csq_constraint_morbid_recessive_distal.tsv \
                                                                             src/visualisation/clinvar_by_csq_constraint.py
	source activate ukb
	python3 -m src.visualisation.clinvar_by_csq_constraint data/statistics/clinvar_by_csq_constraint_morbid_recessive_distal.tsv data/plots/clinvar/clinvar_acmg_by_csq_and_constraint_recessive_distal.png

# Plot PTV ascertainment per region
data/plots/clinvar/clinvar_ptvs_ascertainment.png : data/statistics/clinvar_ptv_ascertainment.tsv \
                                                    src/visualisation/clinvar_ptv_ascertainment.py
	source activate ukb
	python3 -m src.visualisation.clinvar_ptv_ascertainment

# Plot ACMG classifications of PTVs in ClinVar by region and constraint
data/plots/clinvar/clinvar_acmg_by_region_and_constraint.png : data/statistics/clinvar_acmg_by_region_and_constraint.tsv \
                                                               src/visualisation/clinvar_acmg_by_region_and_constraint.py
	source activate ukb
	python3 -m src.visualisation.clinvar_acmg_by_region_and_constraint

# Plot phyloP, AlphaMissense, and pext scores
data/plots/orthogonal_metrics/phylop_alphamis_pext.png : data/statistics/orthogonal_metrics.tsv.gz \
                                                         src/visualisation/orthogonal_metrics.py
	source activate ukb
	python3 -m src.visualisation.orthogonal_metrics

# Plot CADD scores
data/plots/orthogonal_metrics/cadd.png : data/statistics/orthogonal_metrics_cadd.tsv.gz \
                                         src/visualisation/orthogonal_metrics_cadd.py
	source activate ukb
	python3 -m src.visualisation.orthogonal_metrics_cadd

# MAPS
data/plots/maps/maps.png : data/statistics/maps.tsv \
                           src/visualisation/maps.py
	source activate ukb
	python3 -m src.visualisation.maps