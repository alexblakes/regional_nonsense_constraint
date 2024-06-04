#!/usr/bin/env bash
set -euo pipefail

# Concatenate the raw Roulette mutation rate scores
## NB the raw Roulette scores lack a "chr" prefix

export TMP="data/tmp"
export FILE_HEADER="data/manual/header_lines_roulette_raw_bgz.txt"
FILE_OUT="data/interim/roulette/raw_cds_combined.vcf.gz"
FILE_PATHS_FILTERED="${TMP}/roulette_filtered_paths.txt"
FILE_PATHS_RAW="${TMP}/roulette_raw_paths.txt"
export REGIONS="data/interim/gencode_v39_canonical_cds_no_chr.bed"
ROULETTE_DIR="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/roulette/GRCh38/vcfs"

mkdir -p $TMP

# Write the paths of the raw Roulette scores to a text file.
find $ROULETTE_DIR -type f -name *.vcf.bgz | sort -k1,1V > $FILE_PATHS_RAW

# Write the new VCF header to a text file.
cat << EOF > $FILE_HEADER
##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1,length=248956422,assembly=GRCh38>
##contig=<ID=2,length=242193529,assembly=GRCh38>
##contig=<ID=3,length=198295559,assembly=GRCh38>
##contig=<ID=4,length=190214555,assembly=GRCh38>
##contig=<ID=5,length=181538259,assembly=GRCh38>
##contig=<ID=6,length=170805979,assembly=GRCh38>
##contig=<ID=7,length=159345973,assembly=GRCh38>
##contig=<ID=8,length=145138636,assembly=GRCh38>
##contig=<ID=9,length=138394717,assembly=GRCh38>
##contig=<ID=10,length=133797422,assembly=GRCh38>
##contig=<ID=11,length=135086622,assembly=GRCh38>
##contig=<ID=12,length=133275309,assembly=GRCh38>
##contig=<ID=13,length=114364328,assembly=GRCh38>
##contig=<ID=14,length=107043718,assembly=GRCh38>
##contig=<ID=15,length=101991189,assembly=GRCh38>
##contig=<ID=16,length=90338345,assembly=GRCh38>
##contig=<ID=17,length=83257441,assembly=GRCh38>
##contig=<ID=18,length=80373285,assembly=GRCh38>
##contig=<ID=19,length=58617616,assembly=GRCh38>
##contig=<ID=20,length=64444167,assembly=GRCh38>
##contig=<ID=21,length=46709983,assembly=GRCh38>
##contig=<ID=22,length=50818468,assembly=GRCh38>
##INFO=<ID=PN,Number=1,Type=String,Description="Pentanucleotide context">
##INFO=<ID=MR,Number=1,Type=Float,Description="Roulette mutation rate estimate">
##INFO=<ID=AR,Number=1,Type=Float,Description="Adjusted Roulette mutation rate estimate">
##INFO=<ID=MG,Number=1,Type=Float,Description="gnomAD mutation rate estimate (Karczewski et al. 2020)">
##INFO=<ID=MC,Number=1,Type=Float,Description="Carlson mutation rate estimate (Carlson et al. 2018)">
##FILTER=<ID=high,Description="High quality site.">
##FILTER=<ID=low,Description="Low quality regions as determined by gnomAD sequencing metrics. Mappability<0.5;overlap with>50nt simple repeat;ReadPosRankSum>1;0 SNVs in 100bp window.">
##FILTER=<ID=SFS_bump,Description="Pentamer context with abnormal SFS. The fraction of high-frequency SNVS [0.0005<MAF<=0.2] is greater than 1.5x mutation rate controlled average. Tends to be repetitive contexts.">
##FILTER=<ID=TFBS,Description="Transcription factor binding site as determined by overlap with ChIP-seq peaks.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF

limit_to_regions() {
    # Reheader the VCF. Only keep sites in CDS of canonical transcripts.

    local F_OUT="${TMP}/$(basename -s .vcf.bgz ${1}).vcf.gz"
    cat <(cat ${FILE_HEADER}; tabix -R ${REGIONS} ${1}) | bgzip > $F_OUT
}

# Filter in parallel
export -f limit_to_regions
parallel --arg-file $FILE_PATHS_RAW -j80% limit_to_regions {}

# Concatenate the files
find $TMP -type f -name *.vcf.gz | sort -V > $FILE_PATHS_FILTERED
bcftools concat --file-list $FILE_PATHS_FILTERED -o $FILE_OUT -W=tbi

# Clean up
rm -r $TMP