.ONESHELL :
.PHONY : all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

gene_lists = data/final/gene_list_all.txt \
			 data/final/gene_list_gnomad_constrained.txt \
			 data/final/gene_list_nmd_target_constrained.txt \
			 data/final/gene_list_start_proximal_constrained.txt \
			 data/final/gene_list_long_exon_constrained.txt \
             data/final/gene_list_distal_constrained.txt

constrained_gene_lists = data/final/gene_list_gnomad_constrained.txt \
                         data/final/gene_list_nmd_target_constrained.txt \
						 data/final/gene_list_start_proximal_constrained.txt \
						 data/final/gene_list_long_exon_constrained.txt \
						 data/final/gene_list_distal_constrained.txt

all : $(gene_lists) \
      data/final/gene_list_combined_gprofiler.txt \

# Get constrained gene lists
$(gene_lists) : data/raw/gnomad.v4.0.constraint_metrics.tsv \
                data/interim/gene_ids.tsv \
                data/final/regional_nonsense_constraint.tsv \
				src/gene_ontology/gene_lists.py
	python3 -m src.gene_ontology.gene_lists

# Combine gene lists for gProfiler
data/final/gene_list_combined_gprofiler.txt :
	awk 'FNR==1 {print ">" FILENAME;} {print $0}' $(constrained_gene_lists) \
		> data/final/gene_list_combined_gprofiler.txt

# Next