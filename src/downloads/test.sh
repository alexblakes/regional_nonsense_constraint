#!/bin/bash --login

# Download test file

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# touch test.txt
wget -O - cheat.sh/tail | gzip -c > data/raw/tail