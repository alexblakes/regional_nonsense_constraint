# Regional nonsense constraint
A public repository containing analysis workflow for regional nonsense constraint annotations for protein-coding genes.

## Citation
Please cite our preprint:

Blakes, A. J. M., Whiffin, N., Johnson, C. A., Ellingford, J. & Banka, S. [Regional nonsense constraint offers clinical and biological insights into rare genetic disorders](https://doi.org/10.1101/2024.10.10.24315185). 2024.10.10.24315185 (2024).

## Workflow
Our workflow is organised into multiple experiments, each with its own directory under `src`. The order in which these experiments were run is defined in the `Makefile` in the root directory.

The order in which scripts and modules were run for each experiment is defined in the `Makefile` in the directory of that experiment, for example `src/omim/Makefile`.

## Dependencies
Our software dependencies are defined in .yml files under the `env` directory. Dependencies and virtual environments were managed with Conda. For each step of the analysis, the necessary conda environment is defined in the Makefiles described above. If no Conda environment is given, the default environment used was `env/ukb.yml`.

## Data
Key data resources are described in our README at `data/raw/_README.md`.  

Much of our analysis depends on the gnomAD v4.1 sequencing resources, available for [download here](https://gnomad.broadinstitute.org/data).

Our analysis of *de novo* variation was performed in trio sequencing data from the 100,000 Genomes Project (100KGP), the UK Genomic Medicine Service (GMS), and the [Deciphering Developmental Disorders](https://www.ddduk.org/) study.

Data from the 100KGP and the GMS was accessed through the [National Genomic Research Library](https://www.genomicsengland.co.uk/blog/genomics-101-what-is-the-national-genomic-research-library). Access to the NGRL is available to registered academics with approved projects, as described [here](https://www.genomicsengland.co.uk/join-us).

## Key modules
