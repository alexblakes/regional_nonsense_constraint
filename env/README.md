# Environments documentation
This directory contains the necessary data to reproduce the software environment in which the UKB constraint analysis was performed. 

## YAML files
.yml files allow users to reproduce the conda environments in which the analysis was performed.

## rap_clone_git.sh
A utility script to clone this GitHub repo into a UKB RAP worker.

## rap_env_setup.sh
A utility script for one-time set up of the conda environments on a UKB RAP worker. After running in a DNANexus Jupyter Notebook, the environment should be saved as an image for ease of reuse.