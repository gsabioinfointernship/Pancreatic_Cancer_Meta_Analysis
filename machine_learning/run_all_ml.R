# Top-level runner for Machine Learning module

# Step 1: Run per-cancer models (diagnostic and prognostic)
message("STEP 1: Training per-cancer models...")
source("machine_learning/01_per_cancer_ml.R")

# Step 2: Run pan-cancer aggregation
message("\nSTEP 2: Aggregating results pan-cancer...")
source("machine_learning/02_pancancer_ml.R")
