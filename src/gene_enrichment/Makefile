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
                   data/plots/gene_set_enrichment/gnomad_hpo.svg 

all : $(gene_lists) \
      data/final/gene_list_combined_gprofiler.txt \
	  data/statistics/gene_set_enrichment.tsv \
	  $(enrichment_plots)

# Get constrained gene lists
$(gene_lists) : data/raw/gnomad.v4.0.constraint_metrics.tsv \
                data/interim/gene_ids.tsv \
                data/final/regional_nonsense_constraint.tsv \
				src/gene_enrichment/gene_lists.py
	python3 -m src.gene_enrichment.gene_lists

# Combine gene lists for gProfiler web tool
data/final/gene_list_combined_gprofiler.txt :
	awk 'FNR==1 {print ">" FILENAME;} {print $0}' $(constrained_gene_lists) \
		> data/final/gene_list_combined_gprofiler.txt

# Quantify gene set enrichment
data/statistics/gene_set_enrichment.tsv : src/gene_enrichment/gene_set_enrichment.py \
                                          $(gene_lists)
	# Run this from CSF; it requires API access to gProfiler.
	python3 -m src.gene_enrichment.gene_set_enrichment

# Plot gene set enrichment figures
$(enrichment_plots) : src/gene_enrichment/gse_plots.py \
                      data/statistics/gene_set_enrichment.tsv
	python3 -m src.gene_enrichment.gse_plots

# Next