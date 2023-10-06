"""Defines project-wide constants."""

# Directories
## Data
DATA_DIR = "data"
EXTERNAL_DIR = f"{DATA_DIR}/external"
FINAL_DIR = f"{DATA_DIR}/final"
INTERIM_DIR = f"{DATA_DIR}/interim"
LOGS_DIR = f"{DATA_DIR}/logs"
PLOTS_DIR = f"{DATA_DIR}/plots"
RAW_DIR = f"{DATA_DIR}/raw"
STATISTICS_DIR = f"{DATA_DIR}/statistics"

# Files
## Data
### External
CLINVAR_VARIANT_SUMMARY = f"{EXTERNAL_DIR}/variant_summary.txt"
G2P_CARDIAC = f"{EXTERNAL_DIR}/CardiacG2P_23_8_2023.csv"
G2P_DD = f"{EXTERNAL_DIR}/DDG2P_23_8_2023.csv"
G2P_EYE = f"{EXTERNAL_DIR}/EyeG2P_23_8_2023.csv"
G2P_SKELETAL = f"{EXTERNAL_DIR}/SkeletalG2P_23_8_2023.csv"
G2P_SKIN = f"{EXTERNAL_DIR}/SkinG2P_23_8_2023.csv"
GENCODE_GTF = f"{EXTERNAL_DIR}/gencode.v39.annotation.gtf"
GNOMAD_LOEUF_CONSTRAINT = f"{EXTERNAL_DIR}/supplementary_dataset_11_full_constraint_metrics.tsv"
OMIM_GENEMAP = f"{EXTERNAL_DIR}/genemap2.txt"

### Interim
CANONICAL_CDS_BED = f"{INTERIM_DIR}/gencode_v39_canonical_cds.bed"