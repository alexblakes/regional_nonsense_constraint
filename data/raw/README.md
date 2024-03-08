# Data description
Description of the datasets in the ./data/raw/ directory.

## all.baselevel.021620.tsv
Base-level pext scores.
Downloaded with: qsub src/downloads/download_pext.sh
md5sum (not confirmed) d366d60e9e50ff2287ecd418219cb827

## AlphaMissense_hg38.tsv
SNV-level AlphaMissense scores
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) 178b56ce32b5529f0d28a46c03301b20

## CardiacG2P_23_8_2023.csv -->
Text file containing G2P disease gene annotations.
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) 96da2d56677384e71e1c9579aea5ba5e

## DDG2P_23_8_2023.csv
Text file containing G2P disease gene annotations.
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) 817b82ff899c9ba61816a526f63ad6a6

## EyeG2P_23_8_2023.csv
Text file containing G2P disease gene annotations.
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) 7258668829db523a7273e76ad06e4748

## GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
FASTA file for reference human genome. 
Downloaded with: qsub src/downloads/download_reference_fasta.sh
md5sum (not confirmed) a6da8681616c05eb542f1d91606a7b2f

## GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
FASTA index for reference human genome. 
Downloaded with: qsub src/downloads/download_reference_fasta.sh
md5sum (confirmed) 5fddbc109c82980f9436aa5c21a57c61

## gencode.v39.annotation.gtf
GENCODE v39 comprehensive gene annotation (GTF format)
Downloaded with: qsub src/downloads/download_gencode.sh
md5sum (not confirmed) 1efe634b076160b5b9d7027c91a4a518

## genemap2.txt
A tab-delimited file containing OMIM's Synopsis of the Human Gene Map including additional information such as genomic coordinates and inheritance.
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) 16a5c38b1734d48cc1809ca9a2a618a6

## gnomad.v4.0.constraint_metrics.tsv
gnomAD v4 constraint metrics.
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) 612f3899a85a1b7c6513a9dd65b00485

## grch38_cpg_methylation.tsv
Text file containing methylation level annotations for CpG sites.

This data is derived from a publicly available Hail table at: 
- gs://gnomad-nc-constraint-v31-paper/context_prepared.ht

The methylation data was extracted using Hail. This task cannot be done on the University of Manchester computer cluster owing to software availability. Therefore, these scripts were run on a Spark cluster on the UKB Research Analysis Platform, and subsequently downloaded to the UoM CSF. The relevant UKB RAP scripts are here:
- src/downloads/get_ht.sh
- src/downloads/get_methylation.ipynb

md5sum (not confirmed) d1871d99a0b9c3cdb2738c011477c50d

## hg19ToHg38.over.chain.gz
UCSC liftover chain file for hg19 to hg38 annotations.
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) 9267c9fef79b54962da8efadd0ddf6b6

## hg38.cactus241way.phyloP.bw
BigWig file containing base-level 241-way phyloP scores from Zoonomia project.
Downloaded with: qsub src/downloads/download_phylop.sh
md5sum (confirmed) fffa7057e9afa1e177090108775d6418

## mutation_rate_by_context_methyl.txt
Text file with trinucleotide context mutation rates.
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) ed3e599d57dc49a25738a97dbea2c651

## SkeletalG2P_23_8_2023.csv
Text file containing G2P disease gene annotations.
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) 958ed6bbdcb39fa3573362243a69a397

## SkinG2P_23_8_2023.csv
Text file containing G2P disease gene annotations.
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) 01552c7a0681c9dae4ff5227d39eb5ff

## variant_summary.txt
Text file containing ClinVar variant summary data.
Downloaded with: qsub src/downloads/download_small_files.sh
md5sum (not confirmed) 1e9eb39d6f8856809c65b83268e37b56

## whole_genome_SNVs.tsv.gz
CADD scores for all SNVs in GRCh38.
Downloaded with: qsub src/downloads/download_cadd.sh
md5sum (confirmed) 88577a55f1cd519d44e0f415ba248eb9

## whole_genome_SNVs.tsv.gz.tbi
Tabix index for GRCh38 CADD scores.
Downloaded with: qsub src/downloads/download_cadd.sh
md5sum (confirmed) 347df8fac17ea374c4598f4f44c7ce8b