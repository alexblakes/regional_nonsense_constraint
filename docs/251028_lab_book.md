# 28/10/25 Lab book

## Intersecting constrained regions and Pfam domains
After the filtering and deduplication logic, a large number of duplicate enst/region pairs remain.
On manual inspection, these all seem to be regions overlapping identically sized, but differently-named domains. 
E.g. SAM_1/SAM_2, SH3_2/SH3_9.
It will be safe to arbitrarily drop duplicates again, therefore.