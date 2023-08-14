# Add conda channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Update conda
mamba update --force conda -y

# Install packages which directly affect Jupyter Lab
# Also install nb_conda_kernels
mamba install jupyterlab_code_formatter black jupyterlab-lsp python-lsp-server nb_conda_kernels -y
mamba update --all -y

# Install production packages to a new conda environment
# Also install ipykernel
mamba create -n ukb python=3.10 -y
mamba install -n ukb -y tqdm ipykernel scipy scikit-learn gtfparse numpy pandas statsmodels matplotlib seaborn requests dxpy upsetplot
mamba update -n ukb --all -y
