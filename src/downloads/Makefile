.PHONY: small big test all

SHELL = bash

small_files = data/raw/genemap2.txt \
              data/raw/hg19ToHg38.over.chain \
              data/raw/variant_summary.txt \
			  data/raw/CardiacG2P.csv \
              data/raw/SkeletalG2P.csv \
			  data/raw/DDG2P.csv \
			  data/raw/EyeG2P.csv \
      	      data/raw/SkinG2P.csv \
			  data/raw/gnomad.v4.0.constraint_metrics.tsv \
			  data/raw/mutation_rate_by_context_methyl.txt \
			  data/raw/AlphaMissense_hg38.tsv \

big_files = data/raw/hg38.cactus241way.phyloP.bw \
            data/raw/gencode.v39.annotation.gtf \
            data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
	        data/raw/all.baselevel.021620.tsv \
			data/raw/whole_genome_SNVs.tsv.gz \

small : ${small_files}

big : $(big_files)

test : 
	qsub src/downloads/test.sh

all : small big

# Small files
${small_files}:
	qsub src/downloads/download_small_files.sh

# Big files
data/raw/hg38.cactus241way.phyloP.bw : 
	qsub src/downloads/download_phylop.sh

data/raw/gencode.v39.annotation.gtf : 
	qsub src/downloads/download_gencode.sh

data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna : 
	qsub src/downloads/download_reference_fasta.sh

data/raw/all.baselevel.021620.tsv :
	qsub src/downloads/download_pext.sh

data/raw/whole_genome_SNVs.tsv.gz :
	qsub src/downloads/download_cadd.sh