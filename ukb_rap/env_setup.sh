# This script sets up the conda environments on the UKB cloud worker.

# Update Unix environment
apt update
apt upgrade -y
apt install nano -y

# Install github CLI
type -p curl >/dev/null || ( apt update &&  apt install curl -y)
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg |  dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg \
&&  chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
&& echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" |  tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
&&  apt update \
&&  apt install gh -y

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
