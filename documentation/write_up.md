# Formating tricks

> This is a block indent. It is used for table and figure legends.

`
This is a code block.
`

<div class="alert alert-block alert-info">
    This is an <b>alert block</b>. It is used for visible comments.
</div>

A hidden comment is included below.

[//]: # "This is a hidden comment. It is not included in the preview or output."
    
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

### Variant identification
Variants were identified from coding exome sequencing (CES)[^REF] data in the UK Biobank (UKB)[^REF]. The sequencing methods have been previously described[^REF]. Variants were filtered through the gnomAD sample and variant quality control (QC) pipeline[^REF]. Single nucleotide variants (SNVs) in our transcripts of interest, which passed all variant filters (FILTER="PASS") and had a minimum allele count of 1, were extracted from the 421,212 samples which passed sample-level QC. In total, 10,836,767 SNVs meeting these criteria were identified.

<div class="alert alert-block alert-info">
    The gnomAD QC pipeline may require more description. Is it identical to the ExAC pipeline?
</div>

### Modelling the expected number of variants
This analysis was based on the methods described by the gnomAD group [^REF].

Each genomic position in our transcripts of interest was annotated with its trinucleotide context (the reference base on the + strand, plus the bases immediately adjacent to it) using BCFtools **[version]** and reference genome **[reference genome]**. Each possible variant context (64 trinculeotide contexts * 3 possible alternative alleles = 192 variant contexts) was annotated with the per-base, per-generation mutation rate ("mutability") described in Chen et al. 2022 [^REF].

Every possible coding SNV in our transcripts of interest were annotated with the Variant Effect Predictor (VEP) v105 [^REF] (Ensembl v105, consistent with MANE v0.95 and GENCODE v39). For each SNV, the consequence was annotated against all available transcripts. Subsequently, the VEP output was filtered to our transcripts of interest.

<div class="alert alert-block alert-info">
    This is not the case currently. Instead, I have used the VEP --pick, --pick_order, and --coding_only flags. This is a potential source of error. When we refactor the code and rerun the analysis, I should run this as above.
</div>

To model the expected number of variants in any arbitrary sequence in this cohort, we constructed a weighted least-squares regression model describing the relationship between mutability and the proportion of possible synonymous variants observed for each variant context. The model was weighted by the number of possible synonymous variants in each variant context. CpG transitions (accounting for ~3% of all possible synonymous variants) were excluded from this model and from the remainder of the constraint analysis.

<div class="alert alert-block alert-info">
    More discussion on this will be given in the supplementary information. I think a sensible approach would be to build separate models for non-CpG and CpG variants, and to apply them separately to each transcript.
</div>

To calculate the expected number of variants in each transcript and region, we summed the mutability of every possible variant context for each transcript and NMD region (described above). We then applied our model to predict the expected number of synonymous, missense, and nonsense variants per transcript and NMD region. 

### Regional nonsense constraint
To identify constrained transcripts and regions, we performed a one-sided Z-test comparing the proportion of nonsense variants observed (obs) with the proportion of nonsense variants expected (exp) for each transcript and region. We tested the null hypothesis that the proportion of observed nonsense variants is equal to or greater than the proportion of expected nonsense variants in each transcript and region.

Transcripts and regions which are poorly covered may appear depleted of nonsense variants (false positives). Because coverage summary statistics were not available, we used the proportion of synonymous variants observed as a proxy for coverage. We excluded those transcripts and regions in which the proportion of synonymous variants observed was more than one standard deviation below the proportion of synonymous variants expected (synonymous Z score < -1). 

<div class="alert alert-block alert-info">
    The current synonymous filter is for >= 10 synonymous variants and synonymous Z-score > -1. I should amend this to just the Z-score.
</div>

After correcting for multiple testing with the Benjamini-Hochberg (false discovery rate (FDR)) method, we defined constrained transcripts and regions as those with fewer nonsense variants than expected (one-tailed Z test, P < 0.001), or those with 0 nonsense variants observed and P < 0.01. This second group represents transcripts at the limit of our power to detect constraint. To acheive a significance level of P < 0.01 or P < 0.001 (prior to FDR correction), a minimum of ~5 or ~9 expected nonsense variants are required, respectively. To exclude very large transcripts with modest constraint but highly significant P-values, we excluded transcripts and regions where the ratio of observed : expected variants (O/E) was larger than 0.35. 

<div class="alert alert-block alert-info">
    I still need to apply FDR correction, and decide on appropriate significance cutoffs.  
    Should I apply the O/E cutoff prior to FDR correction?
</div>

### Constrained transcripts in ClinVar
<div class="alert alert-block alert-info">
    Write up later. The current approach is a rough draft.
</div>

### Constrained transcripts in PanelApp
<div class="alert alert-block alert-info">
    Write up later. The current approach is a rough draft.
</div>

## Results
### 38.7% of the coding exome is potentially subject to NMD escape
We identified coding positions in which a premature termination codon (PTC) may evade NMD in 19,982 canonical human transcripts. Using four previously-described NMD escape rules [^REF], we annotated positions which were start-proximal (<150nt downstream of the translation start codon), at the 5' end of long exons (>400nt upstream of a splice donor site), subject to the 50nt rule (within the most 3' 50nt of the penultimate exon) or in the final exon of the transcript. Positions in the latter two groups are collectively refered to as "distal NMD escape" positions. In total, 13,214,649 coding positions (38.7% of the coding exome) are subject to one or more of these NMD escape rules ([Table 1](#nmd_table), [Figure 1](#nmd_figure))

<a name="nmd_table"></a>

|NMD region|Number of positions|% of coding genome|
|----------|------------------:|---:|
|NMD target|20,963,322|61.3%|
|Distal|7,518,656|23.2%|
|Start proximal|2,989,495|8.8%|
|Long exon|2,706,498|8.7%|
|**Total**|**34,177,971**|**100%**|

> **Table 1:** The number of coding positions in each NMD region. Note that some positions are duplicated because they overlap with more than one transcript.

<a name="nmd_figure"></a>

<img src="../plots/230327 Transcript diagram.png" style="width: 800px;"/>

> **Figure 1:** Transcript diagram illustrating NMD escape regions. Thick blue boxes represent coding exons. Dark blue dashes depict NMD escape regions. Labels indicate the percentage of coding bases contained within each region, and the number of transcripts which are constrained for nonsense variants in each region.  

### NNN canonical transcripts are highly intolerant to nonsense variants
To examine selective constraint against nonsense variants at the transcript level, we trained a variant expectation model [^REF] on coding exome sequencing data from 421,212 individuals in the UK Biobank [^REF]. The model, based on the mutation rate of SNVs in different trinucleotide contexts, explains **[95.9%]** of the variance in the proportion of observed rare synonymous variants exome-wide [^SF].

We applied this model to predict the number of synonymous, missense, and nonsense variants expected in this cohort in each canonical transcript. We used a one-sided Z test to test the difference between the proportion of observed and expected variants for each transcript and variant consequence (see Methods).

After excluding poorly covered transcripts and correcting for multiple testing, we identified 2,272 transcripts which were significantly constrained for nonsense variants.

<div class="alert alert-block alert-info">
    Should I apply the O/E < 0.35 cutoff at this stage? Or is this more relevant for clinical variant filtering? 
</div>

<a name="constraint_in_transcripts_by_csq"></a>

<img src="../plots/constraint_in_transcripts_by_csq_combined.png" style="width: 800px;"/>

> **Figure 2:** Transcript-level constraint in 421,212 individuals in the UK Biobank. **Top** The number of expected and observed variants in 19,623 canonical transcripts. The grey dashed line represents x=y, with a slope of 1. The solid blue line is the line of best fit (least squares). The reduced number of observed missense and nonsense variants in many transcripts implies negative selection against these variant types. **Middle** The distribution of observed / expected (O/E) variants per transcript, stratified by variant consequence. The grey dashed line marks O/E = 1. Missense variants are moderately skewed left. Nonsense variants are strongly skewed left. A small peak at the extreme left of the synonymous and missense distributions likely represents transcripts which were poorly covered by sequencing. **Bottom** The distribution of constraint Z scores per transcript , stratified by variant consequence. The grey dashed line marks Z = 0. Vertical red lines mark different P value thresholds for a one-sided Z test (prior to FDR correction). A negative Z score indicates that the proportion of variants observed is lower than expected. 


> Figure: Global Z score distribution  
> Figure: Nonsense Z scores by region  
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
[^SF]: Supplementary_figure

