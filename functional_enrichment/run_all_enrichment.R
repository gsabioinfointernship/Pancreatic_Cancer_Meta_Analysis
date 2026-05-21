# Per Cancer Enrichment Analysis
message("Starting functional enrichment pipeline...")
source("functional_enrichment/01_per_cancer_enrichment.R")

source("functional_enrichment/02_pancancer_enrichment.R")
message("\nAll enrichment analyses complete.")
