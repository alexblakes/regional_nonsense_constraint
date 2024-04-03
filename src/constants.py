"""Defines project-wide constants."""

# Plotting and visualisation
## Labels and naming
NMD_REGIONS = ["nmd_target", "start_proximal", "long_exon", "distal"]

NMD_REGION_LABELS = ["NMD target", "Start proximal", "Long exon", "Distal"]

REGIONS = [
    "transcript",
    "nmd_target",
    "start_proximal",
    "long_exon",
    "distal",
]

REGION_LABELS = [
    "Whole transcript",
    "NMD target",
    "Start proximal",
    "Long exon",
    "Distal",
]

MAPS_CONSEQUENCES = [
    "synonymous_variant",
    "missense_variant",
    "stop_gained",
    "nmd_target",
    "start_proximal",
    "long_exon",
    "distal",
][::-1] # Reversed for plotting

MAPS_LABELS = [
    "Synonymous",
    "Missense",
    "Nonsense (Whole transcript)",
    "Nonsense (NMD target)",
    "Nonsense (Start proximal)",
    "Nonsense (Long exon)",
    "Nonsense (Distal)",
][::-1] # Reversed for plotting

## Plotting
CM = 1 / 2.54  # cm to inches conversion
STYLE_DEFAULT = "src/visualisation/styles/default.mplstyle"
COLOR_VIBRANT = "src/visualisation/styles/color/vibrant.mplstyle"
COLOR_REGIONS = "src/visualisation/styles/color/regions_divergent.mplstyle"
COLOR_MAPS = "src/visualisation/styles/color/maps_consequences.mplstyle"

# Directories
## RDS
RDS_DIR = "/mnt/bmh01-rds/Ellingford_gene"
GNOMAD_DIR = f"{RDS_DIR}/public_data_resources/gnomad/v4.0"

## Data
DATA = "data"
FINAL = f"{DATA}/final"
INTERIM = f"{DATA}/interim"
VEP = f"{INTERIM}/vep_all_snvs"
LOGS = f"{DATA}/logs"
PLOTS = f"{DATA}/plots"
RAW = f"{DATA}/raw"
STATS = f"{DATA}/statistics"

# Files
## GEL exports
DNMS = f"{RDS_DIR}/gel_exports/dnms_for_export_2.tsv"

## RDS
GNOMAD_COVERAGE = f"{GNOMAD_DIR}/coverage/gnomad.exomes.v4.0.coverage.summary.tsv"

## Data
### Raw
ALPHA_MISSENSE = f"{RAW}/AlphaMissense_hg38.tsv"
CLINVAR_VARIANT_SUMMARY = f"{RAW}/variant_summary.txt"
G2P_CARDIAC = f"{RAW}/CardiacG2P_23_8_2023.csv"
G2P_DD = f"{RAW}/DDG2P_23_8_2023.csv"
G2P_EYE = f"{RAW}/EyeG2P_23_8_2023.csv"
G2P_SKELETAL = f"{RAW}/SkeletalG2P_23_8_2023.csv"
G2P_SKIN = f"{RAW}/SkinG2P_23_8_2023.csv"
GENCODE_GTF = f"{RAW}/gencode.v39.annotation.gtf"
GNOMAD_LOEUF_CONSTRAINT = (
    f"{RAW}/supplementary_dataset_11_full_constraint_metrics.tsv"
)
GNOMAD_V4_CONSTRAINT = f"{RAW}/gnomad.v4.0.constraint_metrics.tsv"
GNOMAD_NC_MUTABILITY = f"{RAW}/mutation_rate_by_context_methyl.txt"
GNOMAD_NC_METHYLATION = f"{RAW}/grch38_cpg_methylation.tsv"
OMIM_GENEMAP = f"{RAW}/genemap2.txt"
PEXT_RAW = f"{RAW}/all.baselevel.021620.tsv"
PHYLOP_BW = f"{RAW}/hg38.cactus241way.phyloP.bw"

### Interim
ALPHA_MISSENSE_TIDY = f"{INTERIM}/alpha_missense_tidy.tsv"
CADD_ANNOTATED = f"{INTERIM}/cadd_scores_coding_annotated.tsv"
CANONICAL_CDS_BED = f"{INTERIM}/gencode_v39_canonical_cds.bed"
CANONICAL_CDS_FASTA = f"{INTERIM}/gencode_v39_canonical_cds_seq.tsv"
CANONICAL_CDS_GENE_IDS = f"{INTERIM}/gene_ids.tsv"
CANONICAL_CDS_TRANSCRIPT_IDS = f"{INTERIM}/transcript_ids.tsv"
CDS_ALL_SNVS_TRI_CONTEXT = f"{INTERIM}/cds_trinucleotide_contexts.tsv"
CDS_ALL_SNVS_VCF = f"{INTERIM}/cds_all_possible_snvs.vcf"
CDS_COUNTS_AND_COORDS = f"{INTERIM}/cds_counts_and_coords.tsv"
CDS_PHYLOP_PEXT_MISSENSE = f"{INTERIM}/cds_sites_phylop_pext_missense.tsv"
CLINVAR_LOF_ANNOTATED = "data/interim/clinvar_variants_lof_with_nmd_annotation.tsv"
CLINVAR_SELECTED_TSV = f"{INTERIM}/clinvar_variants_selected.tsv"
CLINVAR_SELECTED_VCF = f"{INTERIM}/clinvar_variants_selected.vcf"
CLINVAR_VEP = "data/interim/clinvar_variants_vep.tsv"
CLINVAR_VEP_TIDY = "data/interim/clinvar_variants_vep_tidy.tsv"
DNMS_ANNOTATED = f"{INTERIM}/dnms_annotated.tsv"
GNOMAD_NC_MUTABILITY_TIDY = f"{INTERIM}/mutation_rate_by_context_methyl_tidy.tsv"
GNOMAD_PASS_SNVS = f"{INTERIM}/gnomad_v4_pass_snvs.tsv"
NMD_ANNOTATIONS = f"{INTERIM}/nmd_annotations.tsv"
_OBS_COUNTS_SYN = f"{INTERIM}/observed_variants_counts_synonymous_cov_"
_OBS_COUNTS_REGIONS = f"{INTERIM}/observed_variants_counts_regions_cov_"
OBS_COUNTS_SYN_20 = f"{INTERIM}/observed_variants_counts_synonymous_cov_20.tsv"
OBS_COUNTS_REGIONS_20 = f"{INTERIM}/observed_variants_counts_regions_cov_20.tsv"
OBS_COUNTS_REGIONS_20_CLEAN = (
    f"{INTERIM}/observed_variants_counts_regions_cov_20_clean.tsv"
)
OMIM_GENEMAP_PARSED = f"{INTERIM}/genemap2_parsed.tsv"
OMIM_GENEMAP_SIMPLE = f"{INTERIM}/genemap2_simple.tsv"
PEXT_BED_37 = f"{INTERIM}/pext_37.bed"
PEXT_BED_38 = f"{INTERIM}/pext_38.bed"
PHYLOP_CDS_SCORES = f"{INTERIM}/phylop_cds_sites.tsv"
PS_SYN_CONTEXT = f"{INTERIM}/proportion_singletons_synonymous_by_context.tsv"
PS_REGIONS = f"{INTERIM}/proportion_singletons_by_csq.tsv"
PS_REGIONS_CPG = f"{INTERIM}/proportion_singletons_by_csq_cpg.tsv"
PS_REGIONS_NON_CPG = f"{INTERIM}/proportion_singletons_by_csq_non_cpg.tsv"
VEP_ALL_SNVS = f"{INTERIM}/cds_all_possible_snvs_vep.vcf"
VEP_ALL_SNVS_TIDY = f"{INTERIM}/cds_all_possible_snvs_vep_tidy.tsv"

### Final
ALL_VARIANTS_MERGED_ANNOTATIONS = f"{FINAL}/all_variants_merged_annotations.tsv"
EXPECTED_VARIANTS_ALL_REGIONS = f"{FINAL}/expected_variants_all_regions.tsv"
GENE_LIST_ALL = f"{FINAL}/gene_list_all.txt"
GENE_LIST_GNOMAD_CST = f"{FINAL}/gene_list_gnomad_constrained.txt"
GENE_LIST_NMD_TARGET = f"{FINAL}/gene_list_nmd_target_constrained.txt"
GENE_LIST_START_PROX = f"{FINAL}/gene_list_start_proximal_constrained.txt"
GENE_LIST_LONG_EXON = f"{FINAL}/gene_list_long_exon_constrained.txt"
GENE_LIST_DISTAL = f"{FINAL}/gene_list_distal_constrained.txt"
MAPS = "data/final/maps.tsv"
REGIONAL_CONSTRAINT_STATS = f"{FINAL}/regional_constraint_stats.tsv"
REGIONAL_NONSENSE_CONSTRAINT = f"{FINAL}/regional_nonsense_constraint.tsv"
PHYLOP_PEXT_MISSENSE_STATS = f"{FINAL}/phylop_pext_missense_annotations_stats.tsv"

### Summary statistics
STATS_NMD_FOOTPRINT = "data/statistics/nmd_footprint.tsv"
STATS_CADD_SYN = f"{STATS}/cadd_synonymous.tsv"
STATS_CADD_MIS = f"{STATS}/cadd_missense.tsv"
STATS_CADD_NON = f"{STATS}/cadd_nonsense.tsv"
STATS_CLINVAR_ASCERTAINMENT = "data/statistics/clinvar_ascertainment.tsv"
STATS_CLINVAR_ACMG_REGION = "data/statistics/clinvar_acmg_by_region.tsv"
STATS_CLINVAR_VUS_REGION = "data/statistics/clinvar_vus_by_region.tsv"
STATS_MAPS = "data/statistics/maps.tsv"
STATS_OE = "data/statistics/oe_transcripts.tsv"
STATS_Z_LOEUF = "data/statistics/z_loeuf.tsv"
STATS_Z_REGIONS = "data/statistics/z_scores_per_region.tsv"
STATS_UPSET_CONSTRAINT = "data/statistics/upset_constrained_regions.tsv"
STATS_PHYLOP_MISSENSE_PEXT = "data/statistics/phylop_missense_pext_scores.tsv"
STATS_GENE_SET_ENRICHMENT = "data/statistics/gene_set_enrichment.tsv"

### Logs
VEP_ALL_SNVS_TIDY_LOG = f"{LOGS}/vep_all_snvs_tidy.log"
