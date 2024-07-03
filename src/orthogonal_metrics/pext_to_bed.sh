#!/usr/bin/env bash
set -euo pipefail

LOG="$(basename $0 .sh).log"
FILE_IN=data/raw/all.baselevel.021620.tsv
FILE_OUT=data/interim/pext_37.bed

# Start the log file
echo $(date +%FT%T) > $LOG

# Process the data
#   Drop the header, add a "chr" prefix, reorder the columns
tail -n +2 $FILE_IN \
| awk -F'[\t:]' -v OFS="\t" '
    $NF != "NaN" {
        $1=$1
        $2 = "chr"$2
        END_POS = $3 + 1
        print $2, $3, END_POS, $1, $NF
        }
    ' \
> $FILE_OUT

# Don't drop duplicates; the input is already sorted
# | sort --stable -k1,1V -k2,2n -k4,4V -u --buffer-size=10% --parallel=4 \

# Log the number of entries
echo "Entries in pext data:" >> $LOG
wc -l $FILE_OUT >> $LOG