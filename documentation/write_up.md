# Title
## Materials and Methods
### Transcript set
One representative transcript (the Ensembl_canonical transcript) was selected for every protein-coding gene in the GENCODE v39 annotion [^REF]. In total, 19,982 transcripts were identified. 18,600 of these are the MANE Select transcript (MANE v0.95) [^REF]. The remaining 1,382 transcripts, in genes which lack a MANE annotation, are the APPRIS principal isoform.

### Annotating NMD regions
For each transcript, the following NMD regions were annotated with custom Python scripts:
- Start proximal: The first 150nt downstream of the translation start site
- Long exon: Positions >400nt upstream of an exon-intron junction
- 50nt rule: The most 3' 50nt in the penultimate (second-most 3') exon
- Last exon: Positions in the final (most 3') exon
- NMD target: Any position within a CDS lacking the annotations above  

"Last exon" and "50nt rule" annotations were unified into one category: "Distal". Positions with multiple NMD annotations were assigned one definitive annotation, with this priority: 
start proximal > distal > long exon. 
|NMD region|Number of positions|
|----------|------------------:|
|NMD target|20,963,322|
|Distal|7,518,656|
|Start proximal|2,989,495|
|Long exon|2,706,498|
|Total|34,177,971|
> Note that some of these positions are duplicated across different transcripts.

### Variant identification
Variants were identified from coding exome sequencing (CES)[^REF] data in the UK Biobank (UKB)[^REF]. The sequencing methods have been previously described[^REF]. Variants were filtered through the gnomAD sample and variant quality control (QC) pipeline[^REF]. Single nucleotide variants (SNVs) in our transcripts of interest, which passed all variant filters (FILTER="PASS") and had a minimum allele count of 1, were extracted from the 421,212 samples which passed sample-level QC. In total, 10,836,767 SNVs meeting these criteria were identified.
> The gnomAD QC pipeline may require more description. Is it identical to the ExAC pipeline?

### Modelling the expected number of variants
This analysis was based on the methods described by the gnomAD group [^REF].

Each genomic position in our transcripts of interest was annotated with its trinucleotide context (the reference base on the + strand, plus the bases immediately adjacent to it) using BCFtools **[version]** and reference genome **[reference genome]**. Each possible variant context (64 trinculeotide contexts * 3 possible alternative alleles = 192 variant contexts) was annotated with the per-base, per-generation mutation rate ("mutability") described in Chen et al. 2022 [^REF].

Every possible coding SNV in our transcripts of interest were annotated with the Variant Effect Predictor (VEP) v105 [^REF] (Ensembl v105, consistent with MANE v0.95 and GENCODE v39). For each SNV, the consequence was annotated against all available transcripts. Subsequently, the VEP output was filtered to our transcripts of interest.
> This is not the case currently. Instead, I have used the VEP --pick, --pick_order, and --coding_only flags. This is a potential source of error. When we refactor the code and rerun the analysis, I should run this as above. 

To model the expected number of variants in any arbitrary sequence in this cohort, we constructed a weighted least-squares regression model describing the relationship between mutability and the proportion of possible synonymous variants observed for each variant context. The model was weighted by the number of possible synonymous variants in each variant context. CpG transitions (accounting for ~3% of all possible synonymous variants) were excluded from this model and from the remainder of the constraint analysis.
> See supplementary information for more discussion on this 

To calculate the expected number of variants in each transcript and region, we summed the mutability of every possible variant context for each transcript and NMD region (described above). We then applied our model to predict the expected number of synonymous, missense, and nonsense variants per transcript and NMD region. 

### Regional nonsense constraint
To identify constrained transcripts and regions, we performed a one-sided Z-test comparing the proportion of nonsense variants observed (obs) with the proportion of nonsense variants expected (exp) for each transcript and region. We tested the null hypothesis that the proportion of observed nonsense variants is equal to or greater than the proportion of expected nonsense variants in each transcript and region.

Transcripts and regions which are poorly covered may appear depleted of nonsense variants (false positives). Because coverage summary statistics were not available, we used the proportion of synonymous variants observed as a proxy for coverage. We excluded those transcripts and regions in which the proportion of synonymous variants observed was more than one standard deviation below the proportion of synonymous variants expected. 
> The current synonymous filter is for >= 10 synonymous variants and synonymous Z-score > -1. I should amend this to just the Z-score.

After correcting for multiple testing with the Benjamini-Hochberg (false discovery rate (FDR)) method, we defined constrained transcripts and regions as those with fewer nonsense variants than expected (one-tailed Z-test, P < 0.001), or those with 0 nonsense variants observed and P < 0.01. This second group represents transcripts at the limit of our power to detect constraint. To acheive a significance level of P < 0.01 or P < 0.001 (prior to FDR correction), a minimum of ~5 or ~9 expected nonsense variants are required, respectively. To exclude very large transcripts with modest constraint but highly significant P-values, we excluded transcripts and regions where the ratio of observed : expected variants (O/E) was larger than 0.35. 
> Still need to apply FDR correction, and decide on appropriate significance cutoffs.
> Apply the O/E >0.35 exclusion before FDR correction? This will reduce the severity of the correction. And that is a good thing!

### Constrained transcripts in ClinVar
> Write up later. The current approach is a rough draft.

### Constrained transcripts in PanelApp
> Write up later. The current approach is a rough draft.

## Results
### 38.7% of the coding exome is potentially subject to NMD escape
<img src="../plots/230327 Transcript diagram.png" alt="" style="height: 350px;"/>
> Figure: NMD regions diagram

### NNN canonical transcripts are highly intolerant to nonsense variants

> Figure: Global Z-score distribution
> Supplementary data: Constraint summary statistics  
> Of which, how many have a pLI / LOEUF annotation, and how many are new?  
> - Small genes  

### NNN canonical transcripts exhibit regional nonsense constraint

> Figure: Regional Z-score distributions  
> How many show both global and regional nonsense constraint?
> - NMD is damaging
> - NMD-escape is damaging
> - Both NMD and NMD-escape are damaging


### Limitations
## References
## Supplementary information
### Building the expectation model
### Comparison with other constraint statistics
### Selecting interesting transcript sets

examined the relationship between the mutability of each variant context and the proportion of possible synonymous variants which were observed in the cohort. Synonymous variants were chosen as a class of coding variants which are generally neutral to selection.  

[^REF]: Reference

