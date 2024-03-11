.ONESHELL :
.PHONY : all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

all : data/interim/cadd_scores_coding.tsv \
      

# Filter for CADD scores of coding SNVs
data/interim/cadd_scores_coding.tsv : src/cadd/cadd_coding_snvs.sh
	$(CONDA_ACTIVATE) bio
	bash src/cadd/cadd_coding_snvs.sh
	$(CONDA_ACTIVATE) ukb