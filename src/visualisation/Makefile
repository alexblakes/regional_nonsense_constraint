.ONESHELL:
.PHONY: all

SHELL = bash

all : data/plots/regions/regions_cds_proportions.png \
      data/plots/clinvar/clinvar_vus_by_region.png \
	  data/plots/clinvar/clinvar_ptvs_ascertainment.png \
	  data/plots/constraint/upset_plot.png \

# Plot proportion of CDS per region
data/plots/regions/regions_cds_proportions.png : data/statistics/regions_cds_proportions.tsv \
                                                 src/visualisation/regions_cds_proportions.py
	source activate ukb
	python3 -m src.visualisation.regions_cds_proportions


# Plot proportion VUS per region
data/plots/clinvar/clinvar_vus_by_region.png : data/statistics/clinvar_vus_by_region.tsv \
                                               src/visualisation/clinvar_proportion_vus_by_region.py
	source activate ukb
	python3 -m src.visualisation.clinvar_proportion_vus_by_region

# Plot PTV ascertainment per region
data/plots/clinvar/clinvar_ptvs_ascertainment.png : data/statistics/clinvar_ptv_ascertainment.tsv \
                                                    src/visualisation/clinvar_ptv_ascertainment.py
	source activate ukb
	python3 -m src.visualisation.clinvar_ptv_ascertainment

# Constraint upset plot
data/plots/constraint/upset_plot.png : data/statistics/constraint_upset_data.tsv \
                                       src/visualisation/constraint_upset_plot.py 
	source activate ukb
	python3 -m src.visualisation.constraint_upset_plot