.ONESHELL :
.PHONY : all

SHELL = bash

plots = data/plots/maps/maps.png \
        data/plots/maps/maps.svg \

all : data/interim/maps_snvs_filtered.tsv \
      data/interim/maps.tsv \
	  $(plots)

# Filter SNVs for MAPS
data/interim/maps_snvs_filtered.tsv : data/interim/cds_all_possible_snvs_annotated.tsv.gz \
                                         src/maps/filter_maps_variants.sh
	source activate vep
	bash src/maps/filter_maps_variants.sh
	# Takes about 2 minutes to run

# Calculate MAPS
data/interim/maps.tsv : data/interim/maps_snvs_filtered.tsv \
                        src/maps/maps.py 
	source activate ukb
	python3 -m src.maps.maps

# Plots MAPS
$(plots) &: data/interim/maps.tsv \
            src/maps/plot_maps.py
	source activate ukb
	python3 -m src.maps.plot_maps
