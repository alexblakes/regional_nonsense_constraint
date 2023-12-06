"""Defines project-wide constants."""

# Directories
## RDS
RDS_DIR = "/mnt/bmh01-rds/Ellingford_gene"
GNOMAD_DIR = f"{RDS_DIR}/public_data_resources/gnomad/v4.0"

## Data
DATA_DIR = "data"
FINAL_DIR = f"{DATA_DIR}/final"
INTERIM_DIR = f"{DATA_DIR}/interim"
LOGS_DIR = f"{DATA_DIR}/logs"
PLOTS_DIR = f"{DATA_DIR}/plots"
RAW_DIR = f"{DATA_DIR}/raw"
STATISTICS_DIR = f"{DATA_DIR}/statistics"

# Files
## RDS
GNOMAD_COVERAGE = f"{GNOMAD_DIR}/coverage/gnomad.exomes.v4.0.coverage.summary.tsv"

## Data
### Raw
ALPHA_MISSENSE = f"{RAW_DIR}/AlphaMissense_hg38.tsv"
CLINVAR_VARIANT_SUMMARY = f"{RAW_DIR}/variant_summary.txt"
G2P_CARDIAC = f"{RAW_DIR}/CardiacG2P_23_8_2023.csv"
G2P_DD = f"{RAW_DIR}/DDG2P_23_8_2023.csv"
G2P_EYE = f"{RAW_DIR}/EyeG2P_23_8_2023.csv"
G2P_SKELETAL = f"{RAW_DIR}/SkeletalG2P_23_8_2023.csv"
G2P_SKIN = f"{RAW_DIR}/SkinG2P_23_8_2023.csv"
GENCODE_GTF = f"{RAW_DIR}/gencode.v39.annotation.gtf"
GNOMAD_LOEUF_CONSTRAINT = f"{RAW_DIR}/supplementary_dataset_11_full_constraint_metrics.tsv"
GNOMAD_V4_CONSTRAINT = f"{RAW_DIR}/gnomad.v4.0.constraint_metrics.tsv"
GNOMAD_NC_MUTABILITY = f"{RAW_DIR}/mutation_rate_by_context_methyl.txt"
GNOMAD_NC_METHYLATION = f"{RAW_DIR}/grch38_cpg_methylation.tsv"
OMIM_GENEMAP = f"{RAW_DIR}/genemap2.txt"
PEXT_RAW = f"{RAW_DIR}/all.baselevel.021620.tsv"
PHYLOP_BW = f"{RAW_DIR}/hg38.cactus241way.phyloP.bw"

### Interim
ALPHA_MISSENSE_TIDY = f"{INTERIM_DIR}/alpha_missense_tidy.tsv"
CANONICAL_CDS_BED = f"{INTERIM_DIR}/gencode_v39_canonical_cds.bed"
CANONICAL_CDS_FASTA = f"{INTERIM_DIR}/gencode_v39_canonical_cds_seq.tsv"
CANONICAL_CDS_GENE_IDS = f"{INTERIM_DIR}/gene_ids.tsv"
CDS_ALL_SNVS_TRI_CONTEXT = f"{INTERIM_DIR}/cds_trinucleotide_contexts.tsv"
CDS_ALL_SNVS_VCF = f"{INTERIM_DIR}/cds_all_possible_snvs.vcf"
CDS_COUNTS_AND_COORDS = f"{INTERIM_DIR}/cds_counts_and_coords.tsv"
NMD_ANNOTATIONS = f"{INTERIM_DIR}/nmd_annotations.tsv"
GNOMAD_NC_MUTABILITY_TIDY = f"{INTERIM_DIR}/mutation_rate_by_context_methyl_tidy.tsv"
GNOMAD_PASS_SNVS = f"{INTERIM_DIR}/gnomad_v4_pass_snvs.tsv"
VEP_ALL_SNVS = f"{INTERIM_DIR}/cds_all_possible_snvs_vep.vcf"
VEP_ALL_SNVS_TIDY = f"{INTERIM_DIR}/cds_all_possible_snvs_vep_tidy.tsv"
_OBS_COUNTS_SYN = f"{INTERIM_DIR}/observed_variants_counts_synonymous_cov_"
_OBS_COUNTS_REGIONS = f"{INTERIM_DIR}/observed_variants_counts_regions_cov_"
OBS_COUNTS_SYN_20 = f"{INTERIM_DIR}/observed_variants_counts_synonymous_cov_20.tsv"
OBS_COUNTS_REGIONS_20 = f"{INTERIM_DIR}/observed_variants_counts_regions_cov_20.tsv"
OBS_COUNTS_REGIONS_20_CLEAN = f"{INTERIM_DIR}/observed_variants_counts_regions_cov_20_clean.tsv"
PEXT_BED_37 = f"{INTERIM_DIR}/pext_37.bed"
PEXT_BED_38 = f"{INTERIM_DIR}/pext_38.bed"
PHYLOP_CDS_SCORES = f"{INTERIM_DIR}/phylop_cds_sites.tsv"
PS_SYN_CONTEXT = f"{INTERIM_DIR}/proportion_singletons_synonymous_by_context.tsv"
PS_REGIONS = f"{INTERIM_DIR}/proportion_singletons_by_csq.tsv"
PS_REGIONS_CPG = f"{INTERIM_DIR}/proportion_singletons_by_csq_cpg.tsv"
PS_REGIONS_NON_CPG = f"{INTERIM_DIR}/proportion_singletons_by_csq_non_cpg.tsv"

### Final
ALL_VARIANTS_MERGED_ANNOTATIONS = f"{FINAL_DIR}/all_variants_merged_annotations.tsv"
EXPECTED_VARIANTS_ALL_REGIONS = f"{FINAL_DIR}/expected_variants_all_regions.tsv"
REGIONAL_CONSTRAINT_STATS = f"{FINAL_DIR}/regional_constraint_stats.tsv"
REGIONAL_NONSENSE_CONSTRAINT = f"{FINAL_DIR}/regional_nonsense_constraint.tsv"

### Logs
VEP_ALL_SNVS_TIDY_LOG = f"{LOGS_DIR}/vep_all_snvs_tidy.log"

## Visualisation
DEFAULT_MPL = "src/visualisation/default.mplstyle"
CM = 1/ 2.54 # cm to inches conversion