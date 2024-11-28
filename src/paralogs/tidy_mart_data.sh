#!/usr/bin/env bash
# Tidy the paralogue data from BioMart
set -euo pipefail
source src/utils.sh

FILE_LOG="data/logs/src.paralogs.tidy_mart_data.log"
FILE_PARALOGUES="data/raw/ensembl_paralogs_full.txt.gz"
FILE_GENE_IDS="data/final/canonical_gene_ids.txt"
FILE_OUT="data/interim/ensembl_paralogs_filtered.tsv"

# Filter paralogues
#   awk for autosomal genes only
#   awk for genes in the canonical genes list only
< $FILE_PARALOGUES zcat \
| tail -n +2 \
| tee >(MyLog "Entries: " $(wc -l)) \
| tee >(MyLog "Unique ENSG IDs: " $(cut -f2 | sort -u | wc -l)) \
| awk '$1 >= 1 && $1 <= 22 && $3 >= 1 && $3 <= 22' \
| awk 'NR==FNR {ids[$1]; next} ($2 in ids) && ($4 in ids) {print $0}' $FILE_GENE_IDS - \
| tee >(MyLog "Unique ENSG IDs after filtering canonical ENSGs: " $(cut -f2 | sort -u | wc -l)) \
| grep "within_species_paralog" \
| tee >(MyLog "Unique ENSG IDs after filtering for within species paralogs: " $(cut -f2 | sort -u | wc -l)) \
| cut -f 2,4 \
> $FILE_OUT

sleep 1