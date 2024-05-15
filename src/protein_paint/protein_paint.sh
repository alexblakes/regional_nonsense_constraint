#!/usr/bin/env bash
set -euo pipefail

# Get text files for ProteinPaint

FILE_ARGS="data/manual/protein_paint_args.txt"
CLINVAR="data/interim/clinvar_variants_annotated.vcf"

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

# Compress and index the ClinVar VCF
bgzip -kf $CLINVAR
tabix -f "${CLINVAR}.gz"

# Extract variants from ClinVar and gnomAD
parallel --arg-file $FILE_ARGS --colsep ' ' bash src/protein_paint/clinvar.sh {1} {2} {3} {4}
parallel --arg-file $FILE_ARGS --colsep ' ' bash src/protein_paint/gnomad.sh {1} {2} {3} {4} stop_gained
parallel --arg-file $FILE_ARGS --colsep ' ' bash src/protein_paint/gnomad.sh {1} {2} {3} {4} frameshift_variant

# Extract nonsense variants in gnomAD for the PCDH clusters
for CSQ in "stop_gained" "frameshift_variant"; do
    bash src/protein_paint/gnomad_pcdh.sh PCDHA1 chr5 140786136 141012347 ${CSQ} &
    bash src/protein_paint/gnomad_pcdh.sh PCDHGA1 chr5 141330514 141512975 ${CSQ}
    wait
done

# Tidy up
rm $FILE_ARGS