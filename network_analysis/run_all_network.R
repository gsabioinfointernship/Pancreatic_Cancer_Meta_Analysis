message("Starting network analysis pipeline...")

# Per-cancer PPI networks, hub genes, modules
source("network_analysis/01_per_cancer_network.R")

# Pan-cancer network, shared hub analysis, cross-cancer comparison
source("network_analysis/02_pancancer_network.R")

message("\nAll network analyses complete.")
