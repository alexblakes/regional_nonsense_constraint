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
	   data/interim/observed_variants_counts_regions_cov_20_clean.tsv \
	   data/final/expected_variants_all_regions.tsv \
	   data/final/regional_constraint_stats.tsv \
	   data/final/regional_nonsense_constraint.tsv \
	   data/interim/proportion_singletons_synonymous_by_context.tsv \
	   data/interim/proportion_singletons_by_csq.tsv \

# Files which takes several minutes to create
medium : data/interim/cds_counts_and_coords.tsv \
         data/interim/gene_ids.tsv \
         data/interim/nmd_annotations.tsv \
	     data/interim/cds_all_possible_snvs.vcf \
	     data/interim/cds_trinucleotide_contexts.tsv \
		 data/interim/cds_all_possible_snvs_vep_tidy.tsv \
		 data/interim/gnomad_v4_pass_snvs.tsv \
		 data/interim/observed_variants_counts_regions_cov_20.tsv \
		 data/interim/observed_variants_counts_synonymous_cov_20.tsv \
		 data/interim/pext_37.bed \
		 data/interim/pext_38.bed \

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
data/interim/observed_variants_counts_regions_cov_20.tsv \
data/interim/observed_variants_counts_synonymous_cov_20.tsv : data/final/all_variants_merged_annotations.tsv \
                                                              src/data/observed_variants_counts_and_mutability.py 
	python3 -m src.data.observed_variants_counts_and_mutability --coverage 0 10 20 30

# Clean regional variant count data
data/interim/observed_variants_counts_regions_cov_20_clean.tsv : data/interim/observed_variants_counts_regions_cov_20.tsv \
                                                                 src/data/observed_variants_counts_regions_clean.py
	python3 -m src.data.observed_variants_counts_regions_clean

# Get expected variants for all regions
data/final/expected_variants_all_regions.tsv : data/interim/observed_variants_counts_synonymous_cov_20.tsv \
                                               data/interim/observed_variants_counts_regions_cov_20_clean.tsv \
											   src/constraint/expected_variants.py
	python3 -m src.constraint.expected_variants

# Get nonsense-constrained regions
data/final/regional_constraint_stats.tsv \
data/final/regional_nonsense_constraint.tsv : src/constraint/constraint_statistics.py \
                                              src/constraint/regional_nonsense_constraint.py \
											  data/final/expected_variants_all_regions.tsv \
											  data/raw/gnomad.v4.0.constraint_metrics.tsv
	python3 -m src.constraint.constraint_statistics
	python3 -m src.constraint.regional_nonsense_constraint

# Get proportion of singletons for synonymous variant contexts
data/interim/proportion_singletons_synonymous_by_context.tsv : data/final/all_variants_merged_annotations.tsv \
                                                               src/constraint/proportion_singletons_syn_contexts.py
	python3 -m src.constraint.proportion_singletons_syn_contexts

# Get proportion singletons for all variant consequences
data/interim/proportion_singletons_by_csq.tsv : src/constraint/proportion_singletons_syn_contexts.py \
                                                src/constraint/proportion_singletons_all_csqs.py
	python3 -m src.constraint.proportion_singletons_all_csqs

# Get pext scores (hg19) in bed format
data/interim/pext_37.bed : data/raw/all.baselevel.021620.tsv \
                           src/functional_clinical/pext_tidy.py
	python3 -m src.functional_clinical.pext_tidy

# Liftover pext scores to hg38
data/interim/pext_38.bed : data/interim/pext_37.bed \
                           src/functional_clinical/pext_liftover.sh
	$(CONDA_ACTIVATE) bio
	bash src/functional_clinical/pext_liftover.sh
	$(CONDA_ACTIVATE) ukb

# Annotate CDS sites with phyloP scores
data/interim/phylop_cds_sites.tsv : data/raw/hg38.cactus241way.phyloP.bw \
                                    data/interim/gencode_v39_canonical_cds.bed \
									src/functional_clinical/phylop_extract_scores.py
	python3 -m src.functional_clinical.phylop_extract_scores

# Next 