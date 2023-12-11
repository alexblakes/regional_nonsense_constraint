#! usr/bin/env bash

# Combine VEP-annotated SNVs
#
# Combine the VEP annotation of all possible CDS SNVs, such that:
# 1. The header, starting with '#', is removed
# 2. Only variants with "stop_gained", "missense", or "synonymous" consequences are kept
# 3. The "location" column is split into "chr" and "pos" columns
# 4. Consequences are sanitised. For example "missense_variant,splice_region_variant" 
#    becomes "missense_variant".
# 5. Only variants in the canonical transcript list are kept.
#
# Save the result to a TSV.

cat data/interim/vep_all_snvs/out_*.tsv | \
grep -v '^#' | \
grep -f data/interim/transcript_ids.tsv | \
grep -E 'stop_gained|missense|synonymous' | \
awk -F '[:\t]' -v OFS='\t' '{ $1=$1; print $0 }' | \
awk -F '\t' -v OFS='\t' '{ sub(".*missense.*", "missense_variant", $5); print $0 }' | \
awk -F '\t' -v OFS='\t' '{ sub(".*synonymous.*", "synonymous_variant", $5); print $0 }' | \
awk -F '\t' -v OFS='\t' '{ sub(".*stop_gained.*", "stop_gained", $5); print $0 }' > \
data/interim/cds_all_possible_snvs_vep_combined.tsv