.ONESHELL:
.PHONY: all

SHELL = bash

all : data/interim/genemap2_parsed.tsv \
      data/interim/genemap2_simple.tsv \

# Parse genemap2 data
data/interim/genemap2_parsed.tsv : data/raw/genemap2.txt \
                                   src/omim/omim_parse_genemap.py
	source activate ukb
	python3 -m src.omim.omim_parse_genemap

# Simplify genemap2 annotations
data/interim/genemap2_simple.tsv : data/interim/genemap2_parsed.tsv \
                                   src/omim/omim_simplify_genemap.py
	source activate ukb
	python3 -m src.omim.omim_simplify_genemap