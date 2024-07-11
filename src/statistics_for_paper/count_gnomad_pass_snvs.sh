#!/usr/bin/env bash
set -euo pipefail

# Count autosomal pass SNVs in gnomAD

SNVS="data/interim/gnomad_v4.1_pass_snvs.vcf.gz"

zcat $SNVS | grep -v "^#" | wc -l | xargs echo "Pass SNVs in gnomAD:"