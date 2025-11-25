# Data description
Description of the datasets in the ./data/raw/ directory.

## 250121_biomart_ensembl_refseq_transcript_ids.txt.gz
Text file containing Ensembl and RefSeq transcript IDs  
Note that the Ensembl version 105 database was not available at this time  
BioMart query:
- Ensembl BioMart version 113
- Dataset: Human genes (GRCh38.p14)
- Attributes: Transcript stable ID, RefSeq mRNA ID  
Downloaded from BioMart at this link: http://www.ensembl.org/biomart/martresults/56?file=martquery_0121165148_333.txt.gz  
Download date: 21/01/25  
md5sum 040386eb4f34fb751cc91f5e596c957d  

## all.baselevel.021620.tsv
Base-level pext scores.  
Downloaded with: qsub src/downloads/download_pext.sh  
md5sum (not confirmed) d366d60e9e50ff2287ecd418219cb827  

## AlphaMissense_hg38.tsv
SNV-level AlphaMissense scores  
Downloaded with: qsub src/downloads/download_small_files.sh  
md5sum (not confirmed) 178b56ce32b5529f0d28a46c03301b20  

## CardiacG2P_23_8_2023.csv
Text file containing G2P disease gene annotations.  
Downloaded with: qsub src/downloads/download_small_files.sh  
md5sum (not confirmed) 96da2d56677384e71e1c9579aea5ba5e  

## clinvar_20250115.vcf.gz
The ClinVar VCF (GRCh38), and its associated .tbi index  
Downloaded with this command: `wget -P data/raw https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20250115.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20250115.vcf.gz.tbi`  
Download date: 20/01/25  
md5sum: 6d52b16b4147395119a6a8e6a1834b81 (matches md5sum on the FTP site)  

## DDG2P_23_8_2023.csv
Text file containing G2P disease gene annotations.  
Downloaded with: qsub src/downloads/download_small_files.sh  
md5sum (not confirmed) 817b82ff899c9ba61816a526f63ad6a6  

## ensembl_paralogs.txt.gz
Compressed text file containing paralogous human gene IDs.  
Download date 26/11/24.  
Downloaded from Ensembl BioMart with this command:  
```bash 
# Note the single quotes
wget 'https://urldefense.com/v3/__http://www.ensembl.org/biomart/martresults/56?file=martquery_1127105518_778.txt.gz__;!!PDiH4ENfjr2_Jw!Ct8Y5jEu47TCXSa3XSD067XPOGeXhZPi-wMtCR4UQtqO1UJQBfjKPZGFu16i1cu6kHNnyzDR-zIO5MOSF_d1i_YHnoo_AMs$' -O data/raw/ensembl_paralogs_full.txt.gz
```
md5sum: 6b82161100249df48ed1cde2ebcbd33a  

### BioMart query details:
Dataset:
- Ensembl Genes 113
- Human genes (GRCh38.p14)

Filters:
- [None selected]
Attributes:
- Chromosome/scaffold name
- Gene stable ID
- Human paralogue chromosome/scaffold name
- Human paralogue gene stable ID
- Human paralogue homology type

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

## genebayes_shet_estimates.tsv
Text file with posterior means and 95% confidence intervals for gene-level shet estimates.  
Reference: Zeng et al. Nature Genetics 2024 "Bayesian estimation of gene constraint from an evolutionary model with gene features"  
DOI: 10.1038/s41588-024-01820-9  
Downloaded from: https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-024-01820-9/MediaObjects/41588_2024_1820_MOESM4_ESM.xlsx  
Download date 18/11/24.  
Manually converted worksheet 2 to TSV.

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

## gwas_catalog_all_association_v1.0.2.tsv
Text file containing GWAS catalog data  
Downloaded with wget -O data/raw/gwas_catalog_all_association_v1.0.2.tsv https://www.ebi.ac.uk/gwas/api/search/downloads/alternative  
Download date 28/01/25  
md5sum 4e014dc56fc08e72a9243bf2cec41521

## hg19ToHg38.over.chain.gz
UCSC liftover chain file for hg19 to hg38 annotations.  
Downloaded with: qsub src/downloads/download_small_files.sh  
md5sum (not confirmed) 9267c9fef79b54962da8efadd0ddf6b6  

## hg38.cactus241way.phyloP.bw
BigWig file containing base-level 241-way phyloP scores from Zoonomia project.  
Downloaded with: qsub src/downloads/download_phylop.sh  
md5sum (confirmed) fffa7057e9afa1e177090108775d6418  

## MANE.GRCh38.v0.95.summary.txt.gz
MANE v0.95 summary text file  
Downloaded with this command: wget -P data/raw https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.95/MANE.GRCh38.v0.95.summary.txt.gz  
Download date: 22/01/25  
md5sum: c10e75d12ebcd847f421681eb79b396a

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

## UP000005640_9606_domain.bed
Uniprot domain annotations for human proteins, in BED format  
Downloaded 27/10/2025 with this command:  
`wget -P data/raw https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/UP000005640_9606_domain.bed`

## ucscGenePfam.tsv
"Pfam domains in GENCODE genes" track from UCSC.  
Downloaded in its entirety from `https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=genes&hgta_track=ucscGenePfam&hgta_table=ucscGenePfam`

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
