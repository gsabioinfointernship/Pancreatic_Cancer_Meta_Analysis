# Master runner for Regulatory and Clinical Deep-Dives

message("\n================================================================================")
message("RUNNING SCIENTIFIC DEEP-DIVE MODULES")
message("================================================================================\n")

# Step 1: Master Regulator Analysis
message("STEP 1: Transcription Factor (Master Regulator) Analysis...")
source("regulatory_analysis/01_per_cancer_regulatory.R")
source("regulatory_analysis/02_pancancer_regulatory.R")

# Step 2: Clinicopathological Correlation
message("\nSTEP 2: Clinicopathological Correlation (Stage/Progression)...")
source("clinical_analysis/01_per_cancer_clinical.R")
source("clinical_analysis/02_pancancer_clinical.R")

message("\n================================================================================")
message("SCIENTIFIC DEEP-DIVES COMPLETE")
message("================================================================================\n")
