.ONESHELL:
.PHONY: all

SHELL = bash

all : data/plots/pfam_domains/nmd_regions_overlapping_pfam_domains.png \
      data/plots/pfam_domains/nmd_regions_overlapping_pfam_domains.svg

data/interim/pfam_domains.bed : data/raw/ucscGenePfam.tsv \
                                src/protein_domains/pfam_domains_to_bed.sh
	source activate vep
	bash src/protein_domains/pfam_domains_to_bed.sh

data/interim/pfam_domains_merged.bed : data/interim/pfam_domains.bed \
                                       src/protein_domains/merge_pfam_domains_by_name.sh
	source activate vep
	bash src/protein_domains/merge_pfam_domains_by_name.sh

data/interim/nmd_regions_constraint.bed : data/interim/nmd_regions_pre_bed_constraint.tsv \
                                          src/protein_domains/nmd_regions_constraint_to_bed.sh
	source activate vep
	bash src/protein_domains/nmd_regions_constraint_to_bed.sh

data/interim/nmd_regions_pfam_intersect.tsv : data/interim/pfam_domains_merged.bed \
                                              data/interim/nmd_regions_constraint.bed \
											  src/protein_domains/intersect_nmd_regions_with_pfam_domains.sh
	source activate vep
	bash src/protein_domains/intersect_nmd_regions_with_pfam_domains.sh

data/final/pfam_domains_in_nmd_regions.tsv : data/interim/nmd_regions_pfam_intersect.tsv \
                                             src/protein_domains/count_pfam_domains_in_nmd_regions.py
	source activate ukb
	python3 -m src.protein_domains.count_pfam_domains_in_nmd_regions

data/statistics/proportion_nmd_regions_overlapping_pfam_domains.tsv : data/final/pfam_domains_in_nmd_regions.tsv \
                                                                      src/protein_domains/get_overlap_counts.py
	source activate ukb
	python3 -m src.protein_domains.get_overlap_counts

data/plots/pfam_domains/nmd_regions_overlapping_pfam_domains.png \
data/plots/pfam_domains/nmd_regions_overlapping_pfam_domains.svg &: data/statistics/proportion_nmd_regions_overlapping_pfam_domains.tsv \
                                                                    src/protein_domains/plot_protein_domain_proportions.py
	source activate ukb
	python3 -m src.protein_domains.plot_protein_domain_proportions
