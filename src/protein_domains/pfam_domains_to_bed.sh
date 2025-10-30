#!/usr/bin/env bash
# Convert Pfam UCSC data to BED6 format

set -eu

source src/utils.sh # _log

FILE_IN="data/raw/ucscGenePfam.tsv"
FILE_OUT="data/interim/pfam_domains.bed"

< $FILE_IN tail -n +2 \
| tee >(_log "Input lines: " "$(wc -l)") \
| tee >(_log "Top 10 domains:\n" "$(tail -n+2 | cut -f 5 | sort | uniq -c | sort -k1,1nr | head -n 10)") \
| cut -f 2- \
| bedtools bed12tobed6 \
| tee >(_log "Output lines: " "$(wc -l)") \
> $FILE_OUT

sleep 1