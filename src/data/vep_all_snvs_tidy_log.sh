#! usr/bin/env bash

# Log results for VEP-annotated SNVs

# Constants
LOG_FILE=data/logs/vep_all_snvs_tidy.log
DATA_FILE=data/interim/cds_all_possible_snvs_vep_tidy.tsv

# Log results
echo $(date +'%d-%m-%Y %T') > $LOG_FILE

N_VEP=$(grep -v '^#' data/interim/cds_all_possible_snvs_vep.vcf | wc -l)
echo "Number of VEP-annotated SNVs: $N_VEP" >> $LOG_FILE

N_TIDY=$(wc -l $DATA_FILE | awk '{ print $1 }')
echo "Number of SNVs after tidying: $N_TIDY" >> $LOG_FILE

echo "Variant counts:" >> $LOG_FILE
awk -F '\t' '{ print $5 }' $DATA_FILE | sort | uniq -c | sort -nr >> $LOG_FILE

ENST_COUNT=$(awk -F '\t' '{ print $6 }' $DATA_FILE | sort -u | wc -l)
echo "Number of transcripts: $ENST_COUNT" >> $LOG_FILE