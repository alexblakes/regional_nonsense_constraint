.ONESHELL:
.PHONY: all

SHELL = bash

all : data/plots/conservation/constraint_conservation_pct_deltas.svg \
      data/plots/conservation/constraint_conservation_pct_deltas.png

data/interim/phylop_stats_per_region.tsv : data/interim/cds_sites_phylop_pext_missense.tsv \
                                           data/final/regional_nonsense_constraint_for_excel.tsv \
										   src/conservation/get_phylop_statistics_by_region.py
	source activate ukb
	python3 -m src.conservation.get_phylop_statistics_by_region

data/plots/conservation/constraint_conservation_pct_deltas.svg \
data/plots/conservation/constraint_conservation_pct_deltas.png &: data/interim/phylop_stats_per_region.tsv \
                                                                  src/conservation/plot_percentile_deltas.py
	source activate ukb
	python3 -m src.conservation.plot_percentile_deltas
	