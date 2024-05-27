.ONESHELL:
.PHONY: all

SHELL = /bin/bash

gene_ids = data/interim/gene_ids.tsv \
           data/final/canonical_gene_ids.txt \
           data/final/canonical_transcript_ids.txt \
           data/final/canonical_gene_symbols.txt \

indexed_regions = data/final/nmd_annotations_simple.tsv.gz \
                  data/final/nmd_annotations_simple.tsv.gz.tbi \
                  data/manual/header_lines_nmd_annotation.txt \

all : data/interim/gencode_v39_canonical_cds.bed \
      data/interim/gencode_v39_canonical_cds_no_chr.bed \
	  data/interim/cds_counts_and_coords.tsv \
	  $(gene_ids) \
	  data/interim/nmd_annotations.tsv \
	  $(indexed_regions) \

# Extract canonical CDS from GTF file
data/interim/gencode_v39_canonical_cds.bed \
data/interim/gencode_v39_canonical_cds_no_chr.bed &: data/raw/gencode.v39.annotation.gtf \
                                                    src/regions/canonical_cds.py
	source activate ukb
	python3 -m src.regions.canonical_cds
	# Remove "chr" prefix
	sed 's/^chr//' data/interim/gencode_v39_canonical_cds.bed > data/interim/gencode_v39_canonical_cds_no_chr.bed

# Get CDS counts and coordinates
data/interim/cds_counts_and_coords.tsv : data/raw/gencode.v39.annotation.gtf \
                                         src/regions/cds_counts_and_coords.py
	source activate ukb
	python3 -m src.regions.cds_counts_and_coords

# Get gene and transcript IDs
$(gene_ids) &: data/raw/gencode.v39.annotation.gtf \
              src/regions/canonical_cds.py \
              src/regions/canonical_cds_gene_ids.py \
              src/regions/gene_ids.sh
	source activate ukb
	python3 -m src.regions.canonical_cds_gene_ids

# Annotate NMD regions
data/interim/nmd_annotations.tsv : data/interim/cds_counts_and_coords.tsv \
                                   src/regions/nmd_annotations.py
	source activate ukb
	python3 -m src.regions.nmd_annotations
	# Takes about 15 minutes to run

# Simplify and index NMD region annotations
$(indexed_regions) &: src/regions/index_nmd_annotation.sh \
                      data/interim/nmd_annotations.tsv
	source activate vep
	bash src/regions/index_nmd_annotation.sh
	