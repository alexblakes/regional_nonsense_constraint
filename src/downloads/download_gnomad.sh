#!/bin/bash --login

# Download GENCODE v39 GTF.

#$ -cwd
#$ -e data/logs/csf/
#$ -o data/logs/csf/

URLS=(
    cheat.sh/head
    cheat.sh/tail
)

for URL in ${URLS[@]};
    do
        NAME=$( echo $URL | awk -F "/" '{print $NF}' )
        wget $URL
        md5sum $NAME >> test.md5
    done

# md5sum -c test.md5