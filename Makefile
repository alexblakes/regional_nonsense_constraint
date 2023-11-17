.ONESHELL:
.PHONY: downloads fast medium slow notebooks all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

# Raw data downloads
downloads: 
	make -f data/raw/Makefile all

# Files which take less than a minute to create
fast : data/interim/gencode_v39_canonical_cds.bed \
       data/interim/gencode_v39_canonical_cds_seq.tsv \
	   data/interim/mutation_rate_by_context_methyl_tidy.tsv \

# Files which takes several minutes to create
medium : data/interim/cds_counts_and_coords.tsv \
         data/interim/gene_ids.tsv \
         data/interim/nmd_annotations.tsv \
	     data/interim/cds_all_possible_snvs.vcf \
	     data/interim/cds_trinucleotide_contexts.tsv \
		 data/interim/cds_all_possible_snvs_vep_tidy.tsv \
		 data/interim/gnomad_v4_pass_snvs.tsv \
		 data/final/observed_variants_counts_region.tsv \
		 data/final/observed_variants_counts_synonymous.tsv \

# Files which take hours to create
slow : data/interim/cds_all_possible_snvs_vep.vcf \
       data/final/all_variants_merged_annotations.tsv \

# Notebooks
notebooks :
	# Expectation model choices
	for N in 0 10 20 30 ; do \
		papermill notebooks/01_expectation_model_choices.ipynb \
		notebooks/01_expectation_model_choices.ipynb \
		-p coverage $$N ; \
	done

# All files
all : downloads fast medium slow notebooks


# Extract canonical CDS from GTF file
data/interim/gencode_v39_canonical_cds.bed : data/raw/gencode.v39.annotation.gtf \
                                             src/data/canonical_cds.py
	python3 -m src.data.canonical_cds

# Get CDS counts and coordinates
data/interim/cds_counts_and_coords.tsv : data/raw/gencode.v39.annotation.gtf \
                                         src/data/cds_counts_and_coords.py
	python3 -m src.data.cds_counts_and_coords

# Get gene IDs
data/interim/gene_ids.tsv : data/raw/gencode.v39.annotation.gtf \
                            src/data/canonical_cds.py \
							src/data/canonical_cds_gene_ids.py
	python3 -m src.data.canonical_cds_gene_ids

# Annotate NMD regions
data/interim/nmd_annotations.tsv : data/interim/cds_counts_and_coords.tsv \
                                   src/data/nmd_annotations.py
	python3 -m src.data.nmd_annotations

# Get FASTA sequences for CDS regions
data/interim/gencode_v39_canonical_cds_seq.tsv : data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
                                                 data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai \
												 data/interim/gencode_v39_canonical_cds.bed \
												 src/data/get_fasta.sh
	$(CONDA_ACTIVATE) bio
	bash src/data/get_fasta.sh
	$(CONDA_ACTIVATE) ukb

# Get all possible SNVs and trinucleotide contexts
data/interim/cds_all_possible_snvs.vcf :      data/interim/gencode_v39_canonical_cds_seq.tsv \
                                              src/data/coding_snvs.py
	python3 -m src.data.coding_snvs

data/interim/cds_trinucleotide_contexts.tsv : data/interim/gencode_v39_canonical_cds_seq.tsv \
                                              src/data/coding_snvs.py
	python3 -m src.data.coding_snvs

# Extract SNVs from gnomAD
data/interim/gnomad_v4_pass_snvs.tsv : data/interim/gencode_v39_canonical_cds.bed \
                                       src/data/extract_pass_variants.sh \
									   src/data/batch_extract_pass_variants.sh \
									   src/data/combine_pass_snvs.sh
	bash src/data/batch_extract_pass_variants.sh
	qsub src/data/combine_pass_snvs.sh

# Annotate all possible CDS SNVs with VEP
# This script runs over ~ 10 hours
data/interim/cds_all_possible_snvs_vep.vcf : data/interim/cds_all_possible_snvs.vcf \
                                             src/data/vep_all_snvs.sh
	$(CONDA_ACTIVATE) bio
	bash src/data/vep_all_snvs.sh
	$(CONDA_ACTIVATE) ukb

# Tidy the VEP-annotated SNVs and log the results
data/interim/cds_all_possible_snvs_vep_tidy.tsv : data/interim/cds_all_possible_snvs_vep.vcf \
                                                  src/data/vep_all_snvs_tidy.sh
	bash src/data/vep_all_snvs_tidy.sh
	python3 -m src.data.vep_all_snvs_tidy_log

# Tidy mutability data
data/interim/mutation_rate_by_context_methyl_tidy.tsv : src/data/mutability_data.py
	python3 -m src.data.mutability_data

# Merge annotations for all SNVs
data/final/all_variants_merged_annotations.tsv : data/interim/nmd_annotations.tsv \
                                                 data/interim/cds_all_possible_snvs.vcf \
												 data/interim/cds_trinucleotide_contexts.tsv \
												 data/interim/cds_all_possible_snvs_vep_tidy.tsv \
												 data/interim/gnomad_v4_pass_snvs.tsv \
												 data/interim/mutation_rate_by_context_methyl_tidy.tsv \
												 src/data/observed_variants.py
	python3 -m src.data.observed_variants

# Get variant counts and mutability
data/interim/observed_variants_counts_regions_cov_30.tsv \
data/interim/observed_variants_counts_synonymous_cov_30.tsv : data/final/all_variants_merged_annotations.tsv \
                                                              src/data/observed_variants_counts_and_mutability.py 
	python3 -m src.data.observed_variants_counts_and_mutability -c 0 10 20 30

