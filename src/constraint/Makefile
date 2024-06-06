.ONESHELL :
.PHONY : all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

gene_lists = data/final/gene_list_all.txt \
			 data/final/gene_list_gnomad_constrained.txt \
			 data/final/gene_list_nmd_target_constrained.txt \
			 data/final/gene_list_start_proximal_constrained.txt \
			 data/final/gene_list_long_exon_constrained.txt \
             data/final/gene_list_distal_constrained.txt \
			 
all : data/interim/observed_variants_counts_regions_cov_20.tsv \
      data/final/expected_variants_all_regions.tsv \
      $(gene_lists) \
      data/plots/constraint/venn.png \

# Get variant counts and mutation rates
data/interim/observed_variants_counts_regions_cov_20.tsv : data/interim/cds_all_possible_snvs_annotated.tsv.gz \
                                                           src/constraint/observed_variants_counts_and_mutability.py
	source activate ukb
	python3 -m src.constraint.observed_variants_counts_and_mutability --coverage 20

# Get expected variants for all regions
data/final/expected_variants_all_regions.tsv : data/interim/observed_variants_counts_regions_cov_20.tsv \
                                               src/constraint/expected_variants.py
	source activate ukb
	python3 -m src.constraint.expected_variants

# Quantify constraint statistics


# # Get constrained gene lists
# $(gene_lists) : data/raw/gnomad.v4.0.constraint_metrics.tsv \
#                 data/interim/gene_ids.tsv \
#                 data/final/regional_nonsense_constraint.tsv \
# 				src/constraint/gene_lists.py
# 	python3 -m src.constraint.gene_lists

# # Plot Venn diagrams
# data/plots/constraint/venn.png : $(gene_lists) \
#                                  src/constraint/venn_diagrams.py
# 	python3 -m src.constraint.venn_diagrams

# Next