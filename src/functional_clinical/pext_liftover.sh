# Liftover the pext annotations from hg19 to hg38
# Run this script within the "bio" conda environment

liftOver data/interim/pext_37.bed data/raw/hg19ToHg38.over.chain data/interim/pext_38.bed data/logs/pext_liftover_unmapped.txt
