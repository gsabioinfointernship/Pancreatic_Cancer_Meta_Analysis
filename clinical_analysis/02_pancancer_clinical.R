source("clinical_analysis/00_clinical_functions.R")

# 1. Setup
pancancer_csv <- "outputs/pancancer/clinical_analysis/"
pancancer_fig <- "outputs/pancancer/figures/clinical_analysis/"
dir.create(pancancer_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(pancancer_fig, showWarnings = FALSE, recursive = TRUE)

message("\n=== PAN-CANCER CLINICAL ANALYSIS AGGREGATION ===")

# ── Aggregate Risk-Stage Merged Tables
merged_files <- list.files("outputs", pattern = "_risk_clinical_merged.csv", recursive = TRUE, full.names = TRUE)

if (length(merged_files) > 0) {
  all_data <- lapply(merged_files, function(f) {
    cancer_name <- gsub("_risk_clinical_merged.csv", "", basename(f))
    read.csv(f) |> 
      mutate(
        cancer = cancer_name,
        across(!any_of("risk_score"), as.character)
      )
  }) |> bind_rows()
  
  # Convert risk_score back to numeric if it was somehow converted
  all_data$risk_score <- as.numeric(all_data$risk_score)
  
  write.csv(all_data, file.path(pancancer_csv, "pancancer_all_risk_clinical.csv"), row.names = FALSE)
  
  # ── Global Risk vs Stage Boxplot
  p_global <- ggplot(all_data, aes(x = Stage, y = risk_score, fill = Stage)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    facet_wrap(~ cancer, scales = "free_y") +
    scale_fill_brewer(palette = "Reds") +
    labs(
      title = "Pan-Cancer: Risk Score vs. Pathologic Stage",
      x = "Stage", y = "LASSO-Cox Risk Score"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "gray95"),
          strip.text = element_text(face = "bold"))
          
  ggsave(file.path(pancancer_fig, "pancancer_risk_vs_stage_grid.png"), 
         plot = p_global, width = 10, height = 8, dpi = 600)
}

message("  Pan-cancer clinical aggregation complete.")
