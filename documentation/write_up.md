# Title
## Materials and Methods
### Transcript set
For consistency of annotation, one representative transcript (the Ensembl_canonical transcript) was selected for every protein-coding gene in the GENCODE v39 annotion [^REF]. In total, 19,982 transcripts were identified. 18,600 of these are the MANE_Select transcript (MANE v0.95) [^REF]. The remaining 1,382 transcripts, in genes which lack a MANE annotation, are the APPRIS principal isoform.

### Annotating NMD regions
The following NMD regions were annotated with custom Python scripts:
- Start proximal: The first 150nt positions downstream of the translation start site
- Long exon: Positions >400nt upstream of an exon-intron junction
- 50nt rule: The most 3' 50nt in the penultimate (second-most 3') exon
- Last exon: Positions in the final (most 3') exon of the transcript
- NMD-target: Any position within a CDS, with no other NMD annotation  

"Last exon" and "50nt rule" annotations were unified into one category: "Distal". Positions with multiple NMD annotations were assigned one definitive annotation, with this priority: 
start proximal > distal > long exon. 
|NMD region|Number of positions|
|----------|-------------------:|
|NMD target|20,963,322|
|Distal|7,518,656|
|Start proximal|2,989,495|
|Long exon|2,706,498|
|Total|34,177,971|
> Note that some of these positions are duplicated across different transcripts.

### Variant identification
Variants were identified from coding exome sequencing (CES)[^REF] data in the UK Biobank (UKB)[^REF], processed through the gnomAD sample and variant quality control (QC) pipeline[^REF]. The sequencing methods have been previously described[^REF]. Single nucleotide variants (SNVs) in our transcripts of interest, which passed all variant filters (FILTER="PASS") and had a minimum allele count of 1, were extracted from the 421,212 samples which passed sample-level QC. In total, 10,836,767 SNVs meeting these criteria were identified.
> The gnomAD QC pipeline may require more description. Is it identical to the ExAC pipeline?

### Variant consequence annotation
Variants were annotated with VEP v105 (Ensembl v105, consistent with MANE v0.95 and GENCODE v39). For each SNV, the consequence was annotated against all available transcripts. Subsequently, the VEP output was filtered to our transcripts of interest.
> This is not the case currently. Instead, I have used the VEP --pick, --pick_order, and --coding_only flags. This is a potential source of error. When we refactor the code and rerun the analysis, I will run this as above.

### Expecation model
> See supplementary information for more discussion on this 

## Results
## References
## Supplementary information
[^REF]: Reference
