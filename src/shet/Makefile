.ONESHELL:
.PHONY: all

SHELL = bash

plots = data/plots/shet/venn_shet_constraint.png \
        data/plots/shet/venn_shet_constraint.svg \
		data/plots/shet/venn_shet_omim.png \
		data/plots/shet/venn_shet_omim.svg \

all : data/interim/shet_gnomad_regional_constraint.tsv \
      data/interim/shet_gnomad_regional_omim_ad.tsv \
	  $(plots) \

# Merge regional, gnomAD, and shet constraint annotations
data/interim/shet_gnomad_regional_constraint.tsv : data/final/regional_nonsense_constraint.tsv \
                                                   data/interim/gene_ids.tsv \
												   data/raw/genebayes_shet_estimates.tsv \
												   src/shet/merge_regional_shet_loeuf.py
	source activate ukb
	python3 -m src.shet.merge_regional_shet_loeuf

# Find autosomal dominant morbid genes in OMIM
data/interim/shet_gnomad_regional_omim_ad.tsv : data/interim/shet_gnomad_regional_constraint.tsv \
                                                data/interim/genemap2_simple.tsv \
												src/shet/constrained_in_omim.py
	source activate ukb
	python3 -m src.shet.constrained_in_omim

# Plot Venn diagrams
$(plots) &: data/interim/shet_gnomad_regional_constraint.tsv \
            data/interim/shet_gnomad_regional_omim_ad.tsv \
			src/shet/plot_venn.py
	source activate ukb
	python3 -m src.shet.plot_venn data/interim/shet_gnomad_regional_constraint.tsv shet venn_shet_constraint
	python3 -m src.shet.plot_venn data/interim/shet_gnomad_regional_omim_ad.tsv shet venn_shet_omim