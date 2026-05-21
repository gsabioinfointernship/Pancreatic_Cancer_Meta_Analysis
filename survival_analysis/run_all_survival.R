message("Starting survival analysis pipeline...")
message("Using curatedTCGAData + Cox regression on meta-analysis hub genes\n")

# Per-cancer: download TCGA, median-split KM, univariate Cox
source("survival_analysis/01_per_cancer_survival.R")

# Pan-cancer: consistent prognostic genes, HR heatmap, dot plot
source("survival_analysis/02_pancancer_survival.R")

message("\nAll survival analyses complete.")
