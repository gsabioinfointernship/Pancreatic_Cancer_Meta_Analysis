message("Starting recount3 validation pipeline...")

# Per-cancer: download recount3 TCGA, DESeq2, compare with meta-analysis
source("validation/01_per_cancer_validation.R")

# Pan-cancer: cross-cancer summary, consistently validated genes
source("validation/02_pancancer_validation.R")

message("\nAll validation analyses complete.")
