.ONESHELL:
.PHONY: all \
        count_mane_transcripts \
		count_vep_annotations \
		count_gnomad_pass_snvs \
		count_maps_snvs \

SHELL = bash

all : count_mane_transcripts \
      count_vep_annotations \
	  count_gnomad_pass_snvs \
	  count_maps_snvs \

count_mane_transcripts :
	bash src/statistics_for_paper/count_mane_transcripts.sh

count_vep_annotations :
	bash src/statistics_for_paper/count_vep_annotations.sh

count_gnomad_pass_snvs :
	bash src/statistics_for_paper/count_gnomad_pass_snvs.sh

count_maps_snvs :
	bash src/statistics_for_paper/count_maps_snvs.sh