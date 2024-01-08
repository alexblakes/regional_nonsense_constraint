"""Constants used in plotting and data vis."""

# Labels and naming
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

# Plotting
CM = 1 / 2.54  # cm to inches conversion
STYLE_DEFAULT = "src/visualisation/styles/default.mplstyle"
COLOR_VIBRANT = "src/visualisation/styles/color/vibrant.mplstyle"
COLOR_REGIONS = "src/visualisation/styles/color/regions_divergent.mplstyle"
