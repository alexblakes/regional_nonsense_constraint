"""Defines project-wide constants."""

# Directories
## Data
DATA_DIR = "data"
EXTERNAL_DIR = f"{DATA_DIR}/external"
FINAL_DIR = f"{DATA_DIR}/final"
INTERIM_DIR = f"{DATA_DIR}/interim"
PLOTS_DIR = f"{DATA_DIR}/plots"
RAW_DIR = f"{DATA_DIR}/raw"
STATISTICS_DIR = f"{DATA_DIR}/statistics"

# Files
## Data
### External
CLINVAR_VARIANT_SUMMARY = "data/external/variant_summary.txt"
G2P_CARDIAC = "data/external/CardiacG2P_23_8_2023.csv"
G2P_DD = "data/external/DDG2P_23_8_2023.csv"
G2P_EYE = "data/external/EyeG2P_23_8_2023.csv"
G2P_SKELETAL = "data/external/SkeletalG2P_23_8_2023.csv"
G2P_SKIN = "data/external/SkinG2P_23_8_2023.csv"
GENCODE_GTF = "data/external/gencode.v39.annotation.gtf"
GNOMAD_LOEUF_CONSTRAINT = "data/external/supplementary_dataset_11_full_constraint_metrics.tsv"
OMIM_GENEMAP = "data/external/genemap2.txt"

### Interim
CANONICAL_CDS_BED = "data/interim/gencode_v39_canonical_cds.bed"