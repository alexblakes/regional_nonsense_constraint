#!/usr/bin/env bash
set -euo pipefail

# Get text files for ProteinPaint

FILE_ARGS="data/manual/protein_paint_args.txt"
CLINVAR="data/interim/clinvar_variants_annotated.vcf"
COLOURS="data/final/protein_paint/_colours.txt"

# Add entries to the arguments file
cat << EOF > $FILE_ARGS
KANSL1 chr17 46029916 46193429
BTF3 chr5 73498442 73505667
FBXL17 chr5 107859035 108382098
RSPO3 chr6 127118671 127199481
MNX1 chr7 157004854 157010663
SLC16A9 chr10 59650764 59709850
CBX7 chr22 39130772 39152680
MCMBP chr10 119829440 119872843
FOXF1 chr16 86510527 86515422
PROSER1 chr13 39009865 39038089
IRF2BP2 chr1 234604269 234610178
GATA6 chr18 22169589 22202528
AJAP1 chr1 4654609 4792534
PCDHA1 chr5 140786136 141012347
PCDHGA1 chr5 141330514 141512975
EOF

# Compress and index the ClinVar VCF
bgzip -kf $CLINVAR
tabix -f "${CLINVAR}.gz"

# Extract variants from ClinVar and gnomAD
parallel --arg-file $FILE_ARGS --colsep ' ' bash src/protein_paint/clinvar.sh {1} {2} {3} {4}
parallel --arg-file $FILE_ARGS --colsep ' ' bash src/protein_paint/gnomad.sh {1} {2} {3} {4}

# Extract nonsense variants in gnomAD for the PCDH clusters
bash src/protein_paint/gnomad_pcdh.sh PCDHA1 chr5 140786136 141012347 &
bash src/protein_paint/gnomad_pcdh.sh PCDHGA1 chr5 141330514 141512975

# Create text file for colours
echo "N ; #B2182B" > $COLOURS
echo "F ; #f2d0d5" >> $COLOURS

# Tidy up
rm $FILE_ARGS