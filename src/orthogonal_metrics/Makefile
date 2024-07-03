.ONESHELL:
.PHONY: all

SHELL = bash

all : data/interim/pext_37.bed \
      data/interim/pext_38.bed \
	  data/interim/phylop_cds_sites.tsv \

# Get pext scores (hg19) in bed format
data/interim/pext_37.bed : data/raw/all.baselevel.021620.tsv \
                           src/orthogonal_metrics/pext_to_bed.sh
	source activate vep
	bash src/orthogonal_metrics/pext_to_bed.sh
	# Takes about 10 minutes to run

# Liftover pext scores to hg38
data/interim/pext_38.bed : data/interim/pext_37.bed \
                           src/orthogonal_metrics/pext_liftover.sh
	source activate vep
	bash src/orthogonal_metrics/pext_liftover.sh
	# Takes about 10 minutes to run

# Annotate CDS sites with phyloP scores
data/interim/phylop_cds_sites.tsv : data/raw/hg38.cactus241way.phyloP.bw \
                                    data/interim/gencode_v39_canonical_cds.bed \
									src/orthogonal_metrics/phylop_extract_scores.py
	source activate ukb
	python3 -m src.orthogonal_metrics.phylop_extract_scores

# # Tidy AlphaMissense scores
# data/interim/alpha_missense_tidy.tsv : data/raw/AlphaMissense_hg38.tsv \
#                                        src/functional_clinical/alpha_missense_tidy.py
# 	python3 -m src.functional_clinical.alpha_missense_tidy

# # Merge NMD, phyloP, pext, and AlphaMissense annotations
# data/interim/cds_sites_phylop_pext_missense.tsv \
# data/final/phylop_pext_missense_annotations_stats.tsv : data/interim/alpha_missense_tidy.tsv \
#                                                         data/interim/phylop_cds_sites.tsv \
# 														data/interim/pext_38.bed \
# 														src/functional_clinical/merge_orthogonal_annotations.py
# 	python3 -m src.functional_clinical.merge_orthogonal_annotations
