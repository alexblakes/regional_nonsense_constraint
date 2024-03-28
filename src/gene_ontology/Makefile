.ONESHELL :
.PHONY : all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

gene_lists = data/final/gene_list_*.txt

all : $(gene_lists)

# Get constrained gene lists
$(gene_lists) : data/raw/gnomad.v4.0.constraint_metrics.tsv \
			 data/interim/gene_ids.tsv \
			 data/final/regional_nonsense_constraint.tsv \
			 src/gene_ontology/gene_lists.py
	python3 -m src.gene_ontology.gene_lists

# Next