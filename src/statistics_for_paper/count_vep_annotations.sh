#!/usr/bin/env bash
set -euo pipefail

# Count VEP annotations of all possible SNVs

VEP="data/interim/cds_all_possible_snvs_vep.tsv.gz"


zcat $VEP \
| tee \
    >(cut -f 1 | uniq | xargs echo "Contigs:") \
    >(awk '
        {csq[$6]++} 
        END {print "VEP consequence counts:"; for (c in csq) print "\t", c, csq[c]}
    ') \
> /dev/null