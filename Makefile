.ONESHELL:
SHELL = bash

all : data 

data : data/interim/gencode_v39_canonical_cds.bed

# Extract canonical CDS from GTF file
data/interim/gencode_v39_canonical_cds.bed : data/external/gencode.v39.annotation.gtf \
                                             src/data/canonical_cds.py
	python3 -m src.data.canonical_cds

