"""Defines project-wide constants."""

# Directories
## Data
DATA_DIR = "data"
FINAL_DIR = f"{DATA_DIR}/final"
INTERIM_DIR = f"{DATA_DIR}/interim"
LOGS_DIR = f"{DATA_DIR}/logs"
PLOTS_DIR = f"{DATA_DIR}/plots"
RAW_DIR = f"{DATA_DIR}/raw"
STATISTICS_DIR = f"{DATA_DIR}/statistics"

# Files
## Data
### External
CLINVAR_VARIANT_SUMMARY = f"{RAW_DIR}/variant_summary.txt"
G2P_CARDIAC = f"{RAW_DIR}/CardiacG2P_23_8_2023.csv"
G2P_DD = f"{RAW_DIR}/DDG2P_23_8_2023.csv"
G2P_EYE = f"{RAW_DIR}/EyeG2P_23_8_2023.csv"
G2P_SKELETAL = f"{RAW_DIR}/SkeletalG2P_23_8_2023.csv"
G2P_SKIN = f"{RAW_DIR}/SkinG2P_23_8_2023.csv"
GENCODE_GTF = f"{RAW_DIR}/gencode.v39.annotation.gtf"
GNOMAD_LOEUF_CONSTRAINT = f"{RAW_DIR}/supplementary_dataset_11_full_constraint_metrics.tsv"
OMIM_GENEMAP = f"{RAW_DIR}/genemap2.txt"

### Interim
CANONICAL_CDS_BED = f"{INTERIM_DIR}/gencode_v39_canonical_cds.bed"
CANONICAL_CDS_FASTA = f"{INTERIM_DIR}/gencode_v39_canonical_cds_seq.tsv"
CDS_ALL_SNVS_TRI_CONTEXT = f"{INTERIM_DIR}/cds_trinucleotide_contexts.tsv"
CDS_ALL_SNVS_VCF = f"{INTERIM_DIR}/cds_all_possible_snvs.vcf"