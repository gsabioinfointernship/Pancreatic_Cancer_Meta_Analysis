source("clinical_analysis/00_clinical_functions.R")

for (cancer in CANCERS) {
  
  message("\n=== ", toupper(cancer), " ===")
  
  # 1. Setup
  out_csv <- file.path("outputs", cancer, "clinical_analysis")
  out_fig <- file.path("outputs", cancer, "figures", "clinical_analysis")
  dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
  
  # 2. Load Data
  clinical <- load_tcga_clinical(cancer)
  if (is.null(clinical)) {
    message("    SKIPPING — No clinical metadata found.")
    next
  }
  
  risk <- load_risk_scores(cancer)
  if (is.null(risk)) {
    message("    SKIPPING — No ML risk scores found.")
    next
  }
  
  # 3. Clinical Correlation
  message("    Correlating Risk Scores with Tumor Stage...")
  clin_res <- run_clinical_correlation(risk, clinical)
  
  if (is.null(clin_res)) {
    message("    SKIPPING — Stage data not found or insufficient samples.")
    next
  }
  
  # 4. Export & Plot
  write.csv(clin_res, file.path(out_csv, paste0(cancer, "_risk_clinical_merged.csv")), row.names = FALSE)
  
  save_stage_boxplot(clin_res, cancer, file.path(out_fig, paste0(cancer, "_risk_vs_stage_boxplot.png")))
  
  message("    Done: ", cancer)
}

message("\nPer-cancer clinical correlation complete.")
