# Top-level runner for Immune Microenvironment Analysis
# Step 1: Run per-cancer ssGSEA and correlation with Risk Scores
message("STEP 1: Analyzing per-cancer immune infiltration...")
source("immune_analysis/01_per_cancer_immune.R")

# Step 2: Run pan-cancer aggregation
message("\nSTEP 2: Aggregating immune results pan-cancer...")
source("immune_analysis/02_pancancer_immune.R")
