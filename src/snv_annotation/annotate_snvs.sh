#!/usr/bin/env bash
set -euo pipefail

# Annotate SNVs
# The -a files must be in tsv format in order for bcftools annotate to work in a stream.

VEP_VCF="data/interim/cds_all_possible_snvs_vep.vcf.gz"
VEP_TSV="data/interim/cds_all_possible_snvs_vep.tsv.gz"
NMD_TSV="data/final/nmd_annotations_simple.tsv.gz"
GNOMAD_TSV="data/interim/gnomad_v4.1_pass_snvs.tsv.gz"
COVERAGE_TSV="data/interim/gnomad_v4.1_median_coverage.tsv.gz"
ROULETTE_RAW_TSV="data/interim/roulette/raw_cds_combined.tsv.gz"
ROULETTE_SCALED_TSV="data/interim/roulette/scaled_cds_combined_sorted.tsv.gz"
HEADER_VEP="data/manual/header_lines_vep_consequence.txt"
HEADER_NMD="data/manual/header_lines_nmd_annotation.txt"
HEADER_GNOMAD="data/manual/header_lines_gnomad.txt"
HEADER_COVERAGE="data/manual/header_lines_coverage.txt"
HEADER_ROULETTE_RAW="data/manual/header_lines_roulette_raw_tsv.txt"
HEADER_ROULETTE_SCALED="data/manual/header_lines_roulette_scaled.txt"
FASTA="data/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
LOG="data/logs/annotate_snvs.log"
FILE_OUT="data/interim/cds_all_possible_snvs_annotated.tsv.gz"

log() {
    # Log summary data for annotated SNVs.
    awk '
        $8 == "." {COVNULL++}
        $8 == 0 {COV0++}
        $8 > 0 {COV1++}
        $8 >= 10 {COV10++}
        $8 >= 20 {COV20++}
        $8 >= 30 {COV30++}
        $9 != "." {AC++}
        $13 > 0 {MR++}
        $14 > 0 {MR_final++}
        $15 > 0 {MR_scaled ++}
        $16 > 0 {MG++}
        END {print "Number of records: ", NR}
        END {print "Variants where coverage is null: ", COVNULL}
        END {print "Variants with coverage == 0: ", COV0}
        END {print "Variants with coverage > 0: ", COV1}
        END {print "Variants with coverage >= 10: ", COV10}
        END {print "Variants with coverage >= 20: ", COV20}
        END {print "Variants with coverage >= 30: ", COV30}
        END {print "Observed variants: ", AC}
        END {print "Raw Roulette scores: ", MR}
        END {print "Final Roulette scores: ", MR_final}
        END {print "Scaled Roulette scores: ", MR_scaled}
        END {print "gnomAD scores: ", MG}
    ' \
    - > $LOG
}

# Write header line to output file
printf "#CHROM\tPOS\tREF\tALT\tID\tConsequence\tREGION\tCoverage\tAC\tAN\tAF\tPN\tMR\tMR_final\tMR_scaled\tMG\n" | bgzip -f > $FILE_OUT

# Annotate SNVs
bcftools annotate -a $VEP_TSV -c CHROM,POS,REF,ALT,~ID,Consequence -h $HEADER_VEP -Ou $VEP_VCF \
| bcftools annotate -a $NMD_TSV -c CHROM,POS,~ID,REGION -h $HEADER_NMD -Ou \
| bcftools annotate -a $GNOMAD_TSV -c CHROM,POS,REF,ALT,AC,AN,AF -h $HEADER_GNOMAD -Ou \
| bcftools annotate -a $COVERAGE_TSV -c CHROM,POS,Coverage -h $HEADER_COVERAGE -Ou \
| bcftools annotate -a $ROULETTE_RAW_TSV -c CHROM,POS,REF,ALT,-,PN,MR,-,MG,-,MR_final -h $HEADER_ROULETTE_RAW -Ou \
| bcftools annotate -a $ROULETTE_SCALED_TSV -c CHROM,POS,REF,ALT,-,-,MR_scaled -h $HEADER_ROULETTE_SCALED -Ou \
| bcftools query -f '%CHROM %POS %REF %ALT %ID %Consequence %REGION %Coverage %AC %AN %AF %PN %MR %MR_final %MR_scaled %MG' \
| awk -v OFS="\t" '{$1=$1; print $0}' \
| tee >(log) | bgzip -f >> $FILE_OUT

tabix $FILE_OUT -s1 -b2 -e2 -f
