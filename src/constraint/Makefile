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

all : data/interim/oe_stats_regions_cov_20.tsv \
      data/final/expected_variants_all_regions.tsv \
	  data/final/regional_constraint_stats.tsv \
	  data/final/regional_nonsense_constraint.tsv \
      $(gene_lists) \
      $(venn_plots) \
	  $(pairplots) \

# Get variant counts and O/E statistics
data/interim/oe_stats_regions_cov_20.tsv : data/interim/cds_all_possible_snvs_annotated.tsv.gz \
                                           src/constraint/observed_variants_counts_and_oe_stats.py
	source activate ukb
	python3 -m src.constraint.observed_variants_counts_and_oe_stats --coverage 20
	# Takes about 20 minutes to run

# Quantify constraint statistics
data/final/regional_constraint_stats.tsv : data/interim/oe_stats_regions_cov_20.tsv \
                                           src/constraint/get_fdrs_and_gnomad_constraint.py
	source activate ukb
	python3 -m src.constraint.get_fdrs_and_gnomad_constraint

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
                 src/constraint/plot_venn_diagrams.py
	source activate ukb
	python3 -m src.constraint.plot_venn_diagrams

# Plot scatter pairplots
$(pairplots) &: data/final/regional_nonsense_constraint.tsv \
                data/interim/gene_ids.tsv \
				src/constraint/plot_constraint_scatter.py
	source activate ukb
	python3 -m src.constraint.plot_constraint_scatter