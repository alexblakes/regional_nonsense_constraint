.ONESHELL :
.PHONY : all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

all : data/interim/cadd_scores_coding.tsv \
      data/interim/cadd_scores_coding_annotated.tsv \
	  data/statistics/cadd_nonsense.tsv \
	  data/plots/cadd_violins.png \
	  

# Filter for CADD scores of coding SNVs
data/interim/cadd_scores_coding.tsv : src/cadd/cadd_coding_snvs.sh
	$(CONDA_ACTIVATE) bio
	bash src/cadd/cadd_coding_snvs.sh
	$(CONDA_ACTIVATE) ukb

# Annotate coding sites with CADD scores
data/interim/cadd_scores_coding_annotated.tsv : src/cadd/cadd_scores_annotate.py \
                                                data/interim/cadd_scores_coding.tsv
	python3 -m src.cadd.cadd_scores_annotate

# Calculate CADD score statistics
data/statistics/cadd_nonsense.tsv : src/cadd/cadd_stats.py \
                                    data/interim/cadd_scores_coding_annotated.tsv
	python3 -m src.cadd.cadd_stats

# Plot CADD score violins
data/plots/cadd_violins.png : src/cadd/cadd_violins.py \
                              data/interim/cadd_scores_coding_annotated.tsv
	python3 -m src.cadd.cadd_violins

# Next