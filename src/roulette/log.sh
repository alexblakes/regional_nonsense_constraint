#!/usr/bin/env bash
set -euo pipefail

# Log summary data about the processed Roulette data

LOG="data/logs/roulette_scores.log"
RAW="data/interim/roulette/raw_cds_combined.tsv.gz"
SCALED="data/interim/roulette/scaled_cds_combined_sorted.tsv.gz"

printf "Raw scores\n" > $LOG
awk '
    {FILTER[$5]++};
    $7 > 0 {MR++};
    $8 > 0 {AR++};
    $11 > 0 {MRA++};
    $9 > 0 {MG++};
    $10 > 0 {MC++};
    END {print "Raw number of records: ", NR}; 
    END {print "Raw FILTER value counts: "};
    END {for (f in FILTER) print "\t", f, FILTER[f]};
    END {print "Raw scores: " MR}
    END {print "Adjusted scores: " AR}
    END {print "Final Roulette scores: " MRA}
    END {print "gnomAD scores: " MG}
    END {print "Carlson scores: " MC}
' <(zcat $RAW) >> $LOG

printf "\nScaled scores\n" >> $LOG
awk '
    {FILTER[$5]++};
    $6 > 0 {MR++};
    $7 > 0 {SR++};
    END {print "Raw number of records: ", NR}; 
    END {print "Raw FILTER value counts: "};
    END {for (f in FILTER) print "\t", f, FILTER[f]};
    END {print "Raw scores: " MR}
    END {print "Scaled scores: " SR}
' <(zcat $SCALED) >> $LOG
