# Tidy the HMC data
## Get per-site annotations for hg38, drop duplicates and sites missing from hg38

zcat ../data/hmc_wst_uniquevar_outputforshare_v2.txt.gz |\
cut -f 1,3,11 |\
tail -n +2 |\
awk '{FS="\t"}; $2!=""' |\
sort -u >\
../outputs/hmc_38.tsv
