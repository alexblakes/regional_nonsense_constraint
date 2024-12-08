.ONESHELL :
.PHONY : all

SHELL = /bin/bash

gene_lists = data/final/gene_list_all.txt \
			 data/final/gene_list_gnomad_constrained.txt \
			 data/final/gene_list_nmd_target_constrained.txt \
			 data/final/gene_list_start_proximal_constrained.txt \
			 data/final/gene_list_long_exon_constrained.txt \
             data/final/gene_list_distal_constrained.txt \

constrained_gene_lists = data/final/gene_list_gnomad_constrained.txt \
                         data/final/gene_list_nmd_target_constrained.txt \
						 data/final/gene_list_start_proximal_constrained.txt \
						 data/final/gene_list_long_exon_constrained.txt \
						 data/final/gene_list_distal_constrained.txt \

ora_stats = data/statistics/ora_distal.tsv \
            data/statistics/ora_long_exon.tsv \
            data/statistics/ora_nmd_target.tsv \
            data/statistics/ora_start_proximal.tsv \

enrichment_plots = data/plots/gene_set_enrichment/all_gobp.png \
                   data/plots/gene_set_enrichment/all_gobp.svg \
                   data/plots/gene_set_enrichment/all_gomf.png \
                   data/plots/gene_set_enrichment/all_gomf.svg \
                   data/plots/gene_set_enrichment/all_hpo.png \
                   data/plots/gene_set_enrichment/all_hpo.svg \
                   data/plots/gene_set_enrichment/gnomad_gobp.png \
                   data/plots/gene_set_enrichment/gnomad_gobp.svg \
                   data/plots/gene_set_enrichment/gnomad_gomf.png \
                   data/plots/gene_set_enrichment/gnomad_gomf.svg \
                   data/plots/gene_set_enrichment/gnomad_hpo.png \
                   data/plots/gene_set_enrichment/gnomad_hpo.svg \

all : data/final/gene_list_combined_gprofiler.txt \
	  data/statistics/gene_set_enrichment.tsv \
      $(ora_stats) \
	#   $(enrichment_plots) \

# Combine gene lists for gProfiler web tool
data/final/gene_list_combined_gprofiler.txt : $(constrained_gene_lists)
	awk 'FNR==1 {print ">" FILENAME;} {print $0}' $(constrained_gene_lists) \
		> data/final/gene_list_combined_gprofiler.txt

# ORA with WebGestaltR
$(ora_stats) &: src/gene_enrichment/webgestalt.R
	. /etc/profile
	module load apps/gcc/R/4.3.1
	Rscript src/gene_enrichment/webgestalt.R data/final/gene_list_nmd_target_constrained.txt nmd_target
	Rscript src/gene_enrichment/webgestalt.R data/final/gene_list_start_proximal_constrained.txt start_proximal
	Rscript src/gene_enrichment/webgestalt.R data/final/gene_list_long_exon_constrained.txt long_exon
	Rscript src/gene_enrichment/webgestalt.R data/final/gene_list_distal_constrained.txt distal

# # Plot gene set enrichment figures
# $(enrichment_plots) : src/gene_enrichment/gse_plots.py \
#                       data/statistics/gene_set_enrichment.tsv
# 	source activate ukb
# 	python3 -m src.gene_enrichment.gse_plots
