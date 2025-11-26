# Regional nonsense constraint
A public repository describing the workflow to produce regional nonsense constraint annotations for protein-coding genes.

## Citation
Please cite our preprint:

Blakes, A. J. M., Whiffin, N., Johnson, C. A., Ellingford, J. & Banka, S. [Regional nonsense constraint offers clinical and biological insights into rare genetic disorders](https://doi.org/10.1101/2024.10.10.24315185). 2024.10.10.24315185 (2024).

## Workflow
Our workflow is organised into multiple experiments, each with its own directory under `src`. 

The order in which these experiments were run is defined in the `Makefile` in the root directory.

The order in which scripts and modules were run *within each experiment* is defined in the `Makefile` in the directory of that experiment, for example `src/protein_domains/Makefile`.

## Dependencies
Our software dependencies are defined in .yml files under the `env` directory. Dependencies and virtual environments were managed with Conda. For each step of the analysis, the necessary conda environment is defined in the Makefiles described above. The default environment used was `env/ukb.yml`.

## Data
Key data resources are described in our README at `data/raw/_README.md`.  

Much of our analysis depends on the gnomAD v4.1 sequencing resources, available for [download here](https://gnomad.broadinstitute.org/data).

See also our [companion repository](https://github.com/alexblakes/gel_nmd_dnms/) describing our analysis of _de novo_ variants in the Genomics England protected research environment.

## Key modules
- OE95 statistics are calculated within the `src/constraint/observed_variants_counts_and_oe_stats.py` module
- MAPS statistics are calculated within the `src/constraint/maps.py` module
