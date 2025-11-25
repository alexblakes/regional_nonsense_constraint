.ONESHELL:
.PHONY: all

SHELL = bash

all : data/interim/nmd_regions_pre_bed_constraint.tsv \
      data/interim/nmd_regions_color.bed \

data/interim/nmd_regions_pre_bed.tsv : data/interim/nmd_annotations.tsv \
                                       src/constrained_regions_bed/make_regions_bed.sh
	source activate vep
	bash src/constrained_regions_bed/make_regions_bed.sh

data/interim/nmd_regions_pre_bed_constraint.tsv \
data/interim/nmd_regions_color.bed &: data/interim/nmd_regions_pre_bed.tsv \
                                      data/final/regional_nonsense_constraint.tsv \
									  src/constrained_regions_bed/annotate_with_constraint.py
	source activate ukb
	python3 -m src.constrained_regions_bed.annotate_with_constraint
