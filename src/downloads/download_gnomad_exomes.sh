#!/bin/bash

# Script to download gnomAD VCF and TBI for one chromosome.

DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v4.0/vcf"
TBI_NAME=$( echo $1 | awk -F "/" '{print $NF}' )
VCF_NAME=$( echo $2 | awk -F "/" '{print $NF}' )

wget -c -P ${DIR} $1
wget -c -P ${DIR} $2

md5sum ${DIR}/$TBI_NAME >> ${DIR}/tbi.md5
md5sum ${DIR}/$VCF_NAME >> ${DIR}/vcf.md5

# md5sum -c ${DIR}/tbi.md5 ${DIR}/vcf.md5
# md5sum ${DIR}/$TBI_NAME >> tbi.md5
# md5sum -c ${DIR}/tbi.md5