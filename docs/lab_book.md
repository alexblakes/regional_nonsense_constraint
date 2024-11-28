## 18/11/24
GHR (growth hormone receptor) as an interesting case.
It has a long last exon.
The last exon entirely contains the (presumably critical) growth hormone receptor binding domain.
The last exon is strongly depleted for nonsense variants, but does include a number of frameshift variants.
Some of these frameshift variants occur at reasonably high allele frequencies.

"A head-to-head comparison with GeneBayes rather misses the point."

## 25/11/24
How to include an analysis of constraint in paralogs?
It would be nice to find paralogous genes with variable constraint.
For example, gene A is a paralog of gene B. 
gene A has strong distal constraint, but gene B does not.

I think it will be best to group genes into "families" of paralogs.
Or equivalently, assign an identifier to genes which are paralogs of one another.
It will be most straightforward to do this in Python.

## 26/11/24
I have noticed that some ENSG IDs are present only in the left column (target gene), or in the right column (paralog).

This occurs when the paralog falls on chrX or chrY or on a scaffold (i.e. is not autosomal).

I will therefore re-download the data, and filter for autosomal genes only at an earlier stage.

## 27/11/24
Having re-downloaded the data, and filtered for autosomal genes which are in our proscribed gene list, there are still ENSG IDs which are unique to the "targets" set or the "paralogs" set.

Unique to targets:  
family_id   ensg
108         ENSG00000239620
138         ENSG00000144935
184         ENSG00000239620
342         ENSG00000129187
347         ENSG00000120708
360         ENSG00000157833
431         ENSG00000239620
463         ENSG00000213079
556         ENSG00000145919
648         ENSG00000239620
737         ENSG00000239620
774         ENSG00000141469
979         ENSG00000005436
1071        ENSG00000135916
1143        ENSG00000239620
1155        ENSG00000204856
1215        ENSG00000182223
1304        ENSG00000135932
1526        ENSG00000204381
1856        ENSG00000166949
1895        ENSG00000124019

Note that ENSG00000239620 (PRR20G) is seen many times in this list.
It is in the only duplicated ENSG ID in the list.  
Almost all of these are paralog "families" with two members.  
The only exception is family ID 431, which contains PRR20A to PRR20G

The PRR20 genes have these IDs:
ENSG00000227151
ENSG00000229665
ENSG00000204918
ENSG00000234278
ENSG00000204919
ENSG00000239620

All of the family IDs containing PRR20G (ENSG00000239620) are sets of doubles (PRR20G plus another PRR20 gene), or the full family (PRR20A to PRR20G).

Example case: ENSG00000129187  
This ENSG, and its paralog, are each in only one line of the original ensembl data.

**Decision**  
Because of these discrepancies, I will download the data from BioMart again, using no filters.
Let's see if if that improves things.

For several cases (e.g. ENSG00000144935, ENSG00000129187, ENSG00000144935), this re-download appears to have resolved the issue.

I'll run the scripts and check that it is indeed resolved.

## 28/11/24
This issue is resolved after downloading the BioMart data without any pre-applied filters.

### Tasks
- Find genes which are highly constrained in one paralog, but not another.
- E.g. highest variance in constraint scores (this is influenced by the size of the group)
- The greatest absolute difference in constraint scores between any two members of the group (will only highlight extreme examples, not general properties of the group)
- The mean of the largest difference in constraint scores for every member of the group (less sensitive to outliers. In the case of n=2, this is equivalent to the above)

I think the greatest absolute difference in constraint scores will be the most useful statistic.