.ONESHELL:
.PHONY: all

SHELL = bash

gene_ids = data/interim/gene_ids.tsv \
           data/final/canonical_gene_ids.txt \
           data/final/canonical_transcript_ids.txt \

all : data/interim/gencode_v39_canonical_cds.bed \
      data/interim/gencode_v39_canonical_cds_no_chr.bed \
	  data/interim/cds_counts_and_coords.tsv \
	  $(gene_ids)

# Extract canonical CDS from GTF file
data/interim/gencode_v39_canonical_cds.bed \
data/interim/gencode_v39_canonical_cds_no_chr.bed : data/raw/gencode.v39.annotation.gtf \
                                                    src/regions/canonical_cds.py
	python3 -m src.regions.canonical_cds
	# Remove "chr" prefix
	sed 's/^chr//' data/interim/gencode_v39_canonical_cds.bed > data/interim/gencode_v39_canonical_cds_no_chr.bed

# Get CDS counts and coordinates
data/interim/cds_counts_and_coords.tsv : data/raw/gencode.v39.annotation.gtf \
                                         src/regions/cds_counts_and_coords.py
	python3 -m src.regions.cds_counts_and_coords

# Get gene and transcript IDs
$(gene_ids) : data/raw/gencode.v39.annotation.gtf \
              src/regions/canonical_cds.py \
              src/regions/canonical_cds_gene_ids.py \
              src/regions/gene_ids.sh
	python3 -m src.regions.canonical_cds_gene_ids

# # Annotate NMD regions
# data/interim/nmd_annotations.tsv : data/interim/cds_counts_and_coords.tsv \
#                                    src/data/nmd_annotations.py
# 	python3 -m src.data.nmd_annotations
