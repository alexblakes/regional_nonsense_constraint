#!/usr/bin/env bash
source activate vep # Before set command
set -euo pipefail

# Log some key statistics from the SNV annotation pipeline.

TMP="data/tmp"
LOG="data/logs/all_possible_snvs.log"
POSSIBLE_SNVS="data/interim/cds_all_possible_snvs.vcf.gz"
VEP_SNVS="data/interim/cds_all_possible_snvs_vep.tsv.gz"

mkdir -p $TMP

echo "Number of possible SNVs per chromosome:" > $LOG
zcat $POSSIBLE_SNVS | grep -v "^#" | cut -f1 | uniq -c >> $LOG

echo "Total possible SNVs in autosomes:" >> $LOG
zcat $POSSIBLE_SNVS | grep -v "^#" | grep -E "chr[1-9]+" | wc -l >> $LOG

echo "Total possible distinct / unique SNVs in autosomes:" >> $LOG
zcat $POSSIBLE_SNVS | grep -v "^#" | grep -E "chr[1-9]+" | cut -f1,2,4,5 | sort -S10% --parallel 4 -u | wc -l >> $LOG

echo "VEP annotation was done for autosomes only." >> $LOG
echo "VEP annotation dropped, for example, start codon variants." >> $LOG
echo "Sites may overlap multiple canonical transcripts." >> $LOG

echo "SNVs after VEP annotation:" >> $LOG
zgrep -c ^ $VEP_SNVS >> $LOG

echo "Distinct / unique SNVs after VEP annotation:" >> $LOG
zcat $VEP_SNVS | cut -f1-4 | sort -S10% --parallel 4 -u | wc -l >> $LOG