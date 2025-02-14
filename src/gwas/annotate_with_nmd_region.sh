#!/usr/bin/env bash
# Annotate PTVs from the GWAS catalogue with NMD regions.
set -eu

DIR_TMP="data/tmp/gwas"
FILE_IN="data/interim/gwas_catalog_tidy.tsv.gz"
FILE_NMD="data/final/nmd_annotations_simple.tsv.gz"
FILE_NMD_BED="${DIR_TMP}/regions.bed"
FILE_OUT="data/interim/gwas_ptvs_by_region.tsv"

# Create the temporary directory
[[ -d $DIR_TMP ]] && rm -rf $DIR_TMP
mkdir -p $DIR_TMP

# Get NMD regions for autosomal positions
< $FILE_NMD zcat \
| grep "chr[1-9][0-9]*" \
| awk -v OFS="\t" '{print $1, $2-1, $2, $3";"$4}' \
> $FILE_NMD_BED

# bedtools intersect GWAS PTVs with NMD regions
< $FILE_IN zcat \
| tail -n +2 \
| awk -v FS="\t" -v OFS="\t" \
    ' \
        $5 ~ /stop_gained|frameshift_variant/ \
        {print "chr"$2, $3-1, $3, $8";"$5} \
    ' \
| sort -k 1,1V -k2,2n -k3,4 -u \
| bedtools intersect -wa -wb -a $FILE_NMD_BED -b "stdin" \
| awk -v OFS="\t" \
    '{ \
        split($4, a, ";"); \
        split($8, b, ";"); \
        print $1, $2, $3, b[1], a[1], a[2], b[2] \
    }' \
> $FILE_OUT