# Liftover the pext annotations from hg19 to hg38
# Run this script within the "ucsc" conda environment

# Download the pext bed file for hg19
dx download -f outputs/pext_37.bed -o ../outputs/

conda activate ucsc
liftOver ../outputs/pext_37.bed ../data/hg19ToHg38.over.chain.gz ../outputs/pext_38.bed ../outputs/pext_liftover_unmapped.txt
