#!/usr/bin/env bash
set -euo pipefail

# Get text files for ProteinPaint

FILE_ARGS="data/manual/protein_paint_args.txt"

# Overwrite arguments file with first entry
echo "KANSL1 chr17 46029916 46193429" > $FILE_ARGS

# Add other entries
echo "BTF3 chr5 73498442 73505667" >> $FILE_ARGS
echo "FBXL17 chr5 107859035 108382098" >> $FILE_ARGS
echo "RSPO3 chr6 127118671 127199481" >> $FILE_ARGS
echo "MNX1 chr7 157004854 157010663" >> $FILE_ARGS
echo "SLC16A9 chr10 59650764 59709850" >> $FILE_ARGS
echo "CBX7 chr22 39130772 39152680" >> $FILE_ARGS
echo "MCMBP chr10 119829440 119872843" >> $FILE_ARGS
echo "FOXF1 chr16 86510527 86515422 " >> $FILE_ARGS
echo "PROSER1 chr13 39009865 39038089" >> $FILE_ARGS

# # Add genes in protocadherin alpha cluster
# echo "PCDHA1 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA2 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA3 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA4 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA5 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA6 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA7 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA8 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA9 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA10 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA11 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA12 chr5 140786136 141012347" >> $FILE_ARGS
# echo "PCDHA13 chr5 140786136 141012347" >> $FILE_ARGS

while read -r SYMBOL CHROM START END; do
    bash src/protein_paint/clinvar.sh $SYMBOL $CHROM $START $END &
    bash src/protein_paint/gnomad.sh $SYMBOL $CHROM $START $END &
done < $FILE_ARGS

# Combine variants in the PCDHA cluster