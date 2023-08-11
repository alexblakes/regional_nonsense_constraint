# Install packages to a UKB cloud worker using conda.

# There are numerous difficulties managing conda environments on the UKB cloud worker when running the jupyter lab app.
# The best approach is to directly install required packages on to the base conda environment on the UKB worker.
# Do not initiate conda or activate a conda environment before installing these packages.

# Optionally, you can save a snapshot of the environment after these packages are installed.

# I will update this script as new packages are required in UKB.

# First force an update of conda itself
conda update --force conda -y

# Then install desired packages
conda install -c bioconda gtfparse -y
conda install -c anaconda statsmodels scikit-learn -y
conda install -c conda-forge upsetplot -y

# Then update all packages
conda update --all -y
