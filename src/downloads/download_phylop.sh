#!/bin/bash --login

# Download phyloP scores.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

rsync -avz --progress \
        rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus241way/hg38.cactus241way.phyloP.bw \
        data/raw/