#!/bin/bash --login

# Download test file

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

# touch test.txt
wget -P data/raw/ cheat.sh/tail
gzip data/raw/tail