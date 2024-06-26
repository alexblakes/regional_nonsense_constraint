.ONESHELL :
.PHONY : all

SHELL = bash

gene_lists = data/final/gene_list_all.txt \
			 data/final/gene_list_gnomad_constrained.txt \
			 data/final/gene_list_nmd_target_constrained.txt \
			 data/final/gene_list_start_proximal_constrained.txt \
			 data/final/gene_list_long_exon_constrained.txt \
             data/final/gene_list_distal_constrained.txt \

venn_plots = data/plots/constraint/venn.png \
             data/plots/constraint/venn.svg \

pairplots = data/plots/constraint/region_pair_plots.png \
            data/plots/constraint/region_pair_plots.svg \

all : data/interim/observed_variants_counts_regions_cov_20.tsv \
      data/final/expected_variants_all_regions.tsv \
	  data/final/regional_constraint_stats.tsv \
	  data/final/regional_nonsense_constraint.tsv \
      $(gene_lists) \
      $(venn_plots) \
	  $(pairplots) \

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
data/final/regional_constraint_stats.tsv : data/final/expected_variants_all_regions.tsv \
                                           src/constraint/constraint_statistics.py
	source activate ukb
	python3 -m src.constraint.constraint_statistics
	# Takes about 20 minutes to run

# Quantify regional nonsense constraint
data/final/regional_nonsense_constraint.tsv : data/final/regional_constraint_stats.tsv \
                                              src/constraint/regional_nonsense_constraint.py
	source activate ukb
	python3 -m src.constraint.regional_nonsense_constraint

# Get constrained gene lists
$(gene_lists) &: data/final/regional_nonsense_constraint.tsv \
                data/interim/gene_ids.tsv \
				src/constraint/gene_lists.py
	source activate ukb
	python3 -m src.constraint.gene_lists

# Plot Venn diagrams
$(venn_plots) &: $(gene_lists) \
                src/constraint/venn_diagrams.py
	source activate ukb
	python3 -m src.constraint.venn_diagrams

# Plot scatter pairplots
$(pairplots) &: data/final/regional_nonsense_constraint.tsv \
                data/interim/gene_ids.tsv \
				src/constraint/constraint_scatter.py
	source activate ukb
	python3 -m src.constraint.constraint_scatter