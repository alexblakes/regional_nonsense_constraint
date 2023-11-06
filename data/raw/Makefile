.PHONY: big small all test

SHELL = bash

small : genemap2.txt hg19ToHg38.over.chain variant_summary.txt CardiacG2P.csv \
        SkeletalG2P.csv DDG2P.csv EyeG2P.csv SkinG2P.csv

big : hg38.cactus241way.phyloP.bw GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
      GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai gencode.v39.annotation.gtf \
      all.baselevel.021620.tsv

all : small big

test: hg38.cactus241way.phyloP.bw

# Small files


# Big files
hg38.cactus241way.phyloP.bw : 
	qsub src/downloads/download_phylop.sh

GCA_000001405.15_GRCh38_no_alt_analysis_set.fna : 
	qsub src/downloads/download_reference_fasta.sh

gencode.v39.annotation.gtf : 
	qsub src/downloads/download_gencode.sh

all.baselevel.021620.tsv.bgz :
	qsub src/downloads/download_pext.sh