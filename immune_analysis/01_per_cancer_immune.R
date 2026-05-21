source("immune_analysis/00_immune_functions.R")

for (cancer in CANCERS) {
  
  message("\n=== ", toupper(cancer), " ===")
  
  # 1. Setup
  out_csv <- file.path("outputs", cancer, "immune_analysis")
  out_fig <- file.path("outputs", cancer, "figures", "immune_analysis")
  dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
  
  # 2. Load Data
  expr <- load_tcga_expr(cancer)
  if (is.null(expr)) {
    message("    SKIPPING — No TCGA expression data found.")
    next
  }
  
  risk <- load_risk_scores(cancer)
  if (is.null(risk)) {
    message("    SKIPPING — No ML risk scores found.")
    next
  }
  
  # 3. ssGSEA
  ssgsea_res <- run_immune_ssgsea(expr)
  if (is.null(ssgsea_res)) next
  
  # Save raw ssGSEA scores
  saveRDS(ssgsea_res, file.path(out_csv, paste0(cancer, "_immune_ssgsea_scores.rds")))
  
  # 4. Correlation Analysis
  cor_res <- correlate_immune_risk(ssgsea_res, risk)
  if (is.null(cor_res)) {
    message("    SKIPPING — Failed to match samples for correlation.")
    next
  }
  
  # 5. Export & Plots
  write.csv(cor_res$cor_table, file.path(out_csv, paste0(cancer, "_immune_risk_correlation.csv")), row.names = FALSE)
  
  save_immune_risk_barplot(cor_res$cor_table, cancer, 
                            file.path(out_fig, paste0(cancer, "_immune_risk_correlation.png")))
                            
  save_immune_heatmap(ssgsea_res, risk, cancer, 
                       file.path(out_fig, paste0(cancer, "_immune_infiltration_heatmap.png")))
                       
  message("    Done: ", cancer)
}

message("\nPer-cancer immune analysis complete.")
