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
