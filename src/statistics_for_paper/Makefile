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
	  data/statistics/clinvar_ptvs_in_morbid_ad_genes_acmg_by_constraint_and_csq.tsv

count_mane_transcripts :
	bash src/statistics_for_paper/count_mane_transcripts.sh

count_vep_annotations :
	bash src/statistics_for_paper/count_vep_annotations.sh

count_gnomad_pass_snvs :
	bash src/statistics_for_paper/count_gnomad_pass_snvs.sh

count_maps_snvs :
	bash src/statistics_for_paper/count_maps_snvs.sh

data/statistics/clinvar_ptvs_in_morbid_ad_genes_acmg_by_constraint_and_csq.tsv : data/interim/clinvar_variants_constraint_in_morbid_dominant_genes.tsv.gz \
                                                                                 src/statistics_for_paper/count_acmg_by_constrained_region.py
	source activate ukb
	python3 -m src.statistics_for_paper.count_acmg_by_constrained_region
