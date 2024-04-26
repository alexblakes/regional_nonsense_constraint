#!/usr/bin/env bash
set -euo pipefail
# Get list of included genes and transcripts

IN_FILE="data/interim/gene_ids.tsv"
OUT_GENES="data/final/canonical_gene_ids.txt"
OUT_TRANSCRIPTS="data/final/canonical_transcript_ids.txt"

tail -n +2 $IN_FILE | cut -f 1 | sort -u > $OUT_GENES
tail -n +2 $IN_FILE | cut -f 2 | sort -u > $OUT_TRANSCRIPTS

echo "Unique ENSG IDs: $(wc -l < $OUT_GENES)"
echo "Unique ENST IDs: $(wc -l < $OUT_TRANSCRIPTS)"