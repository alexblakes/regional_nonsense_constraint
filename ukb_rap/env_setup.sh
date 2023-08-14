# This script sets up the conda environments on the UKB cloud worker.

# Update conda
conda update --force conda -y

# Install packages which target Jupyter Lab itself
# Including nb_conda_kernels, to allow kernel switching in Jupyter Lab
# Rather than adding channels to the .condarc file, it is quicker to specify
# which channels to install from on the command line.
conda install -c anaconda -y python-lsp-server
conda install -c conda-forge jupyterlab_code_formatter black jupyterlab-lsp nb_conda_kernels -y
conda update --all -y

# Install packages which are used in production into a new conda environment
# Including ipykernel, which allows the new environment to be recognised in Jupyert Lab
conda create -n ukb python=3.8 -y
conda install -n ukb -c bioconda -y dxpy gtfparse
conda install -n ukb -c conda-forge -y tqdm upsetplot 
conda install -n ukb -c anaconda -y ipykernel scipy scikit-learn numpy pandas statsmodels matplotlib seaborn requests
conda update -n ukb --all -y
