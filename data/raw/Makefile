.PHONY: big small all test

SHELL = bash

small_files = data/raw/genemap2.txt data/raw/hg19ToHg38.over.chain \
              data/raw/variant_summary.txt data/raw/CardiacG2P.csv \
              data/raw/SkeletalG2P.csv data/raw/DDG2P.csv data/raw/EyeG2P.csv \
      	      data/raw/SkinG2P.csv

big_files = data/raw/hg38.cactus241way.phyloP.bw data/raw/gencode.v39.annotation.gtf \
            data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
	        data/raw/all.baselevel.021620.tsv

test: qsub src/downloads/test.sh

small:
	qsub src/downloads/download_small_files.sh

big: $(big_files)

all : small big

# Big files
data/raw/hg38.cactus241way.phyloP.bw : 
	qsub src/downloads/download_phylop.sh

data/raw/gencode.v39.annotation.gtf : 
	qsub src/downloads/download_gencode.sh

data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna : 
	qsub src/downloads/download_reference_fasta.sh

data/raw/all.baselevel.021620.tsv :
	qsub src/downloads/download_pext.sh