# ClinVar variants
# This script extracts the relevant data from the ClinVar VCF.
# It writes the output to .tsv.

# Run this line from the "ukb" conda environment
dx download -f -o ../data/ data/clinvar_20230410.vcf.gz data/clinvar_20230410.vcf.gz.tbi 

# Run this line from the "bio" condo environment
conda activate bio
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%CLNREVSTAT\t%CLNSIG\t%CLNVC\t%GENEINFO\t%MC\t%ORIGIN\n' \
    -o ../outputs/clinvar_variants.tsv \
    ../data/clinvar_20230410.vcf.gz

# NB the VCF header annotation for allele origin: Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other