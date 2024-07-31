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
		omim \
		protein_paint \
		statistics_for_plots \
		visualisation \
		figures \
		all \

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

omim :
	make -f src/omim/Makefile all

protein_paint :
	make -f src/protein_paint/Makefile all

statistics_for_plots :
	make -f src/statistics_for_plots/Makefile all

visualisation :
	make -f src/visualisation/Makefile all

figures :

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
	  omim \
	  protein_paint \
	  statistics_for_plots \
	  visualisation \
	  figures \
