#! usr/bin/env bash

# Tidy VEP-annotated SNVs
#
# Tidy the VEP annotation of all possible CDS SNVs, such that:
# 1. The header, starting with '#', is removed
# 2. Only variants with "stop_gained", "missense", or "synonymous" consequences are kept
# 3. The "location" column is split into "chr" and "pos" columns
# 4. Consequences are sanitised. For example "missense_variant,splice_region_variant" 
#    becomes "missense_variant".
#
# Save the result to a TSV.

grep -v '^#' data/interim/cds_all_possible_snvs_vep.vcf | \
grep -E 'stop_gained|missense|synonymous' | \
awk -F '[:\t]' -v OFS='\t' '{ $1=$1; print $0 }' | \
awk -F '\t' -v OFS='\t' '{ sub(".*missense.*", "missense_variant", $5); print $0 }' | \
awk -F '\t' -v OFS='\t' '{ sub(".*synonymous.*", "synonymous_variant", $5); print $0 }' | \
awk -F '\t' -v OFS='\t' '{ sub(".*stop_gained.*", "stop_gained", $5); print $0 }' > \
data/interim/cds_all_possible_snvs_vep_tidy.tsv