#!/bin/bash

# Script to download gnomAD VCF and TBI for one chromosome.

DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v4.1/exomes/vcf"
TBI_NAME=$(basename $1)
VCF_NAME=$(basename $2)

wget -c -P ${DIR} $1
wget -c -P ${DIR} $2

md5sum "${DIR}/${TBI_NAME}" 
md5sum "${DIR}/${VCF_NAME}"