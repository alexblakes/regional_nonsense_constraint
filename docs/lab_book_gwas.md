# Lab book for GWAS analysis

## 03/02/2025 Tidying GWAS catalogue

The unfiltered GWAS catalogue data has 700,076 entries.
After filtering by various criteria, we reduce this to 402,974 entries.
This comprises 243,135 unique SNPs across 24,910 unique traits.
Among these variants, there are only 1,356 nonsense variants (317 unique) and 106 frameshift variants (53 unique)

### Exploring pleiotropy
The idea: for each unique nonsense variant, is the number of impacted traits greater than for synonymous variants for example.
- stop_gained: 317 variants --> 921 traits
- missense_variant: 6,405 variants --> 7,890 traits
- synonymous_variant: 1,787 variants --> 1,571 traits
- intron_variant: 110,721 variants --> 17,213 traits
- intergenic_variant: 105,531 variants --> 13,924 traits

I have explore this more fully with a short script at: src/gwas/explore_pleiotropy.py