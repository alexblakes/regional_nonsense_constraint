#!/usr/bin/env bash
# Merge overlapping Pfam domains with the same name.
# Some Pfam domains with the same name overlap. This causes issues with 
# downstream region annotation. They should be merged.

set -eu

# shellcheck source=/dev/null
source src/utils.sh # _log

FILE_IN="data/interim/pfam_domains.bed"
FILE_OUT="data/interim/pfam_domains_merged.bed"

export dir_tmp
dir_tmp=$(mktemp -d)
trap '[[ -d $dir_tmp ]] && rm -rf $dir_tmp' EXIT

# Split lines by domain name
< $FILE_IN tee >(_log "Input lines: " "$(wc -l)") \
| awk -v dir="$dir_tmp" -v FS="\t" -v OFS="\t" '{print $0 > dir"/"$4".bed"}'

# Merge bed files in parallel
function merge_beds(){
    local file_in=$1
    local bname
    bname=$(basename "$file_in" .bed)
    local file_out="${dir_tmp}/${bname}_merged.bed"
    
    # The merge command keeps domain name and strand information.
    # The awk command restores BED6 format.
    < "$file_in" sort -k1,1V -k2,2n -k3,3n \
    | bedtools merge -s -c 4,6 -o distinct,distinct \
    | awk -v OFS="\t" '$4 = $4OFS"."' \
    > "$file_out"
}
export -f merge_beds

find "$dir_tmp" -type f -name "*.bed" -print0 | parallel -0 merge_beds

# Cat the merged files
find "$dir_tmp" -type f -name "*_merged.bed" -print0 \
| xargs -0 cat \
| tee >(_log "Output lines: " "$(wc -l)") \
> $FILE_OUT 

sleep 1