# Data description
Description of the datasets in the ./data/ directory.

## 230411_biomart_105_transcripts.tsv
.tsv file containing HGNC symbols, ENSG numbers and ENST numbers for all Ensembl_canonical transcripts.
Downloaded from Ensembl BioMart (Ensembl archiver version 105) on 11/04/23.
From this URL: http://dec2021.archive.ensembl.org/biomart/martview/
With these criteria:
- Dataset: 
    - Human genes (GRCh38.p13)
- Filters:
    - Chromosome/scaffold 1-22, X, Y
    - Ensembl Canonical: Only
- Attributes:
    - HGNC symbol
    - Gene stable ID
    - Transcript stable ID
    
## 241-mammalian-2020v2.bigWig
BigWig file containing phyloP scores from Zoonomia project, using multiple alignment of 241 mammals.
This is the updated file (v2).
From this URL: https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/Homo_sapiens/241-mammalian-2020v2.bigWig
Downloaded on 14/06/23.
Note that due to its size (20Gb) this file is kept in the UKB RAP.

## hmc_wst_uniquevar_outputforshare_v2.txt.gz
Gzipped text file containing homologous missense constraint annotations for hg19 and hg38 sites.
Reference: Zhang, X. et al. Genetic constraint at single amino acid resolution improves missense variant prioritisation and gene discovery. medRxiv 2022.02.16.22271023 (2022) doi:10.1101/2022.02.16.22271023.
Downloaded 19/06/23.
From this URL: https://doi.org/10.5281/zenodo.6392153
