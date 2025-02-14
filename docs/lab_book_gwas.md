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

## 04/02/2025 PTVs in GWAS catalogue
For the 402,974 entries in the tidied catalogue, here are some of the most common traits:
                                                   count  rank
trait                                                         
Height                                             16959   1.0
Body mass index                                     7375   2.0
Educational attainment                              4958   3.0
Insomnia                                            3153   7.0
Type 2 diabetes                                     2772  10.0
Schizophrenia                                       2165  17.0
Inguinal hernia                                     2115  19.0
Prostate cancer                                     1259  44.0

There are only five clinical diagnostic phenotypes in the top 50 entries, accounting for ~10,000 associations.

For PTVs (nonsense and frameshift variants) there are only three clinical diagnostic phenotypes in the top 50 traits, spanning 26 associations:
                                                    count  rank
trait                                                          
Body mass index                                        17   1.0
Height                                                 14   2.0
Type 2 diabetes                                        12   3.0
Prostate cancer                                         5  27.0
Type 2 diabetes (adjusted for BMI)                      5  27.0
Age-related hearing impairment                          4  39.0

### Annotation of PTVs in the GWAS catalog
It appears that the consequence annotations in the GWAS catalogue may be inaccurate.
For example, the 19-19268740-C-T allele, which is widely implicated in common liver phenotypes, is annotated as a stop_gained variant. 
It is in fact missense in the canonical transcript of TM6SF2.

According to the GWAS Catalogue documentation, this annotation reflects the strongest consequences of the variant in any transcript.

**NB** Because reference alleles are not given, it is not possible to accurately re-annotate frameshift variants with VEP, for example.

- [ ] Which are the 17 PTVs which are GWAS hits and found in constrained regions?