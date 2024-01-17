# Environments documentation
This directory contains the necessary data to reproduce the software environment in which the UKB constraint analysis was performed. 

## YAML files
.yml files allow users to reproduce the conda environments in which the analysis was performed.

## rap_clone_git.sh
A utility script to clone this GitHub repo into a UKB RAP worker.

## rap_env_setup.sh
A utility script for one-time set up of the conda environments on a UKB RAP worker. After running in a DNANexus Jupyter Notebook, the environment should be saved as an image for ease of reuse.

## pip installs
UpSetPlot 0.9.0 was installed with the following command:
    python3 -m ~/miniforge3/envs/ukb/bin/pip install -U --no-deps -t /mnt/iusers01/bk01/m40482ab/miniforge3/envs/ukb/lib/python3.8/site-packages upsetplot

This command ensures that:
    - Paths are consistent between the conda-installed version of Python and pip
    - No dependencies of UpSetPlot (e.g. Pandas, Matplotlib) are updated
    - Upsetplot is installed to the same directory as packages installed with Mamba