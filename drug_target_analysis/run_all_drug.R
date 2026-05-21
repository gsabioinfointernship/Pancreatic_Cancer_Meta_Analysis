message("Starting drug target analysis pipeline...")

# Per-cancer: DGIdb query, hub gene prioritisation, druggable family classification
source("drug_target_analysis/01_per_cancer_drug.R")

# Pan-cancer: consistent DEG targets, cross-cancer comparison, shared druggable genes
source("drug_target_analysis/02_pancancer_drug.R")

message("\nAll drug target analyses complete.")
