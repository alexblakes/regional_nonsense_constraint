.PHONY: g2p big small all

SHELL = bash

g2p : CardiacG2P.csv SkeletalG2P.csv DDG2P.csv EyeG2P.csv SkinG2P.csv

big : 241-mammalian-2020v2.bigWig GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
      gencode.v39.annotation.gtf all.baselevel.021620.tsv

small : GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai genemap2.txt \
        hg19ToHg38.over.chain variant_summary.txt

all : g2p big small

CardiacG2P.csv : 
	wget -c https://www.ebi.ac.uk/gene2phenotype/downloads/CardiacG2P.csv.gz
	gunzip CardiacG2P.csv.gz

SkeletalG2P.csv :
	wget -c https://www.ebi.ac.uk/gene2phenotype/downloads/SkeletalG2P.csv.gz
	gunzip SkeletalG2P.csv.gz

DDG2P.csv :
	wget -c https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz
	gunzip DDG2P.csv.gz

EyeG2P.csv :
	wget -c https://www.ebi.ac.uk/gene2phenotype/downloads/EyeG2P.csv.gz
	gunzip EyeG2P.csv.gz

SkinG2P.csv :
	wget -c https://www.ebi.ac.uk/gene2phenotype/downloads/SkinG2P.csv.gz
	gunzip SkinG2P.csv.gz

241-mammalian-2020v2.bigWig :
	wget https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/Homo_sapiens/241-mammalian-2020v2.bigWig

GCA_000001405.15_GRCh38_no_alt_analysis_set.fna :
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai :
	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

gencode.v39.annotation.gtf : 
	wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
	gunzip gencode.v39.annotation.gtf.gz

genemap2.txt : 
	wget https://data.omim.org/downloads/dW5qBh3GSkCt5K3BqUQL1w/genemap2.txt

hg19ToHg38.over.chain :
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
	gunzip hg19ToHg38.over.chain.gz

all.baselevel.021620.tsv :
	wget https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-tx-annotation/gnomad_browser/all.baselevel.021620.tsv.bgz -O all.baselevel.021620.tsv.gz
	gunzip all.baselevel.021620.tsv.gz

variant_summary.txt :
	wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
	gunzip variant_summary.txt.gz