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
		notebooks \
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

statistics : 
	papermill src/statistics_for_plots/clinvar_ascertainment.ipynb src/statistics_for_plots/clinvar_ascertainment.ipynb
	papermill src/statistics_for_plots/maps.ipynb src/statistics_for_plots/maps.ipynb
	papermill src/statistics_for_plots/oe.ipynb src/statistics_for_plots/oe.ipynb
	papermill src/statistics_for_plots/z_loeuf.ipynb src/statistics_for_plots/z_loeuf.ipynb
	papermill src/statistics_for_plots/z_distributions.ipynb src/statistics_for_plots/z_distributions.ipynb
	papermill src/statistics_for_plots/upset.ipynb src/statistics_for_plots/upset.ipynb

figures :

notebooks :
	# Expectation model choices
	for N in 0 10 20 30 ; do \
		papermill notebooks/01_expectation_model_choices.ipynb \
		notebooks/01_expectation_model_choices.ipynb \
		-p coverage $$N ; \
	done

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
	  statistics \
	  figures \

