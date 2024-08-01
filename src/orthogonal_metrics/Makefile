.ONESHELL:
.PHONY: all

SHELL = bash

all : data/interim/pext_37.bed \
      data/interim/pext_38.bed \
	  data/interim/phylop_cds_sites.tsv \
	  data/interim/alpha_missense_tidy.tsv \
	  data/interim/cds_sites_phylop_pext_missense.tsv \
	  data/final/phylop_pext_missense_annotations_stats.tsv \
	  data/interim/cadd_scores_coding.tsv \
	  data/interim/cadd_scores_coding_annotated.tsv \

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
	# Takes about 10 minutes to run

# Tidy AlphaMissense scores
data/interim/alpha_missense_tidy.tsv : data/raw/AlphaMissense_hg38.tsv \
                                       src/orthogonal_metrics/alpha_missense_tidy.py
	source activate ukb
	python3 -m src.orthogonal_metrics.alpha_missense_tidy
	# Takes about 10 minutes to run

# Merge NMD, phyloP, pext, and AlphaMissense annotations
data/interim/cds_sites_phylop_pext_missense.tsv \
data/final/phylop_pext_missense_annotations_stats.tsv : data/interim/alpha_missense_tidy.tsv \
                                                        data/interim/phylop_cds_sites.tsv \
														data/interim/pext_38.bed \
														data/interim/nmd_annotations.tsv \
														data/interim/gene_ids.tsv \
														data/final/regional_nonsense_constraint.tsv \
														src/orthogonal_metrics/merge_orthogonal_annotations.py
	source activate ukb
	python3 -m src.orthogonal_metrics.merge_orthogonal_annotations

# Filter for CADD scores of coding SNVs
data/interim/cadd_scores_coding.tsv : src/orthogonal_metrics/cadd_coding_snvs.sh
	source activate vep
	bash src/orthogonal_metrics/cadd_coding_snvs.sh

# Annotate coding sites with CADD scores
data/interim/cadd_scores_coding_annotated.tsv : data/interim/cadd_scores_coding.tsv \
                                                data/interim/nmd_annotations.tsv \
                                                data/final/regional_nonsense_constraint.tsv \
                                                data/interim/cadd_scores_coding.tsv \
												src/orthogonal_metrics/cadd_scores_annotate.py
	source activate ukb
	python3 -m src.orthogonal_metrics.cadd_scores_annotate