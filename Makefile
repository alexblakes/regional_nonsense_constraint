.ONESHELL:
.PHONY: downloads fast medium slow all

SHELL = bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

# Raw data downloads
downloads: 
	make -f data/raw/Makefile all

# Files which take less than a minute to create
fast : data/interim/gencode_v39_canonical_cds.bed \
       data/interim/gencode_v39_canonical_cds_seq.tsv \

# Files which takes several minutes to create
medium : data/interim/cds_counts_and_coords.tsv \
         data/interim/nmd_annotations.tsv \
	     data/interim/cds_all_possible_snvs.vcf \
	     data/interim/cds_trinucleotide_contexts.tsv \
		 data/interim/cds_all_possible_snvs_vep_tidy.tsv \

# Files which take hours to create
slow : data/interim/cds_all_possible_snvs_vep.vcf \

# All files
all : downloads fast medium slow


# Extract canonical CDS from GTF file
data/interim/gencode_v39_canonical_cds.bed : data/raw/gencode.v39.annotation.gtf \
                                             src/data/canonical_cds.py
	python3 -m src.data.canonical_cds

# Get CDS counts and coordinates
data/interim/cds_counts_and_coords.tsv : data/raw/gencode.v39.annotation.gtf \
                                         src/data/cds_counts_and_coords.py
	python3 -m src.data.cds_counts_and_coords

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
	bash src/data/vep_all_snvs_tidy_log.sh
