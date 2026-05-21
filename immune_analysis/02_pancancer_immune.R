source("immune_analysis/00_immune_functions.R")

# 1. Setup pan-cancer output directory
pancancer_csv <- "outputs/pancancer/immune_analysis/"
pancancer_fig <- "outputs/pancancer/figures/immune_analysis/"
dir.create(pancancer_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(pancancer_fig, showWarnings = FALSE, recursive = TRUE)

message("\n=== PAN-CANCER IMMUNE ANALYSIS AGGREGATION ===")

# ── Aggregate Correlation Tables
cor_files <- list.files("outputs", pattern = "_immune_risk_correlation.csv", 
                        recursive = TRUE, full.names = TRUE)

if (length(cor_files) > 0) {
  all_cor <- lapply(cor_files, function(f) {
    cancer_name <- gsub("_immune_risk_correlation.csv", "", basename(f))
    read.csv(f) |> mutate(cancer = cancer_name)
  }) |> bind_rows()
  
  write.csv(all_cor, file.path(pancancer_csv, "pancancer_immune_risk_correlation.csv"), row.names = FALSE)
  
  # Identify universal immune markers of risk
  summary_cor <- all_cor |>
    group_by(Cell_Type) |>
    summarise(
      Mean_r = mean(Spearman_r, na.rm = TRUE),
      SD_r = sd(Spearman_r, na.rm = TRUE),
      N_Sig = sum(P_Value < 0.05, na.rm = TRUE),
      Cancers = paste(cancer[P_Value < 0.05], collapse = "; "),
      .groups = "drop"
    ) |>
    arrange(desc(abs(Mean_r)))
    
  write.csv(summary_cor, file.path(pancancer_csv, "pancancer_consistent_immune_markers.csv"), row.names = FALSE)
  
  # ── Heatmap of Correlations across Cancers
  plot_df <- all_cor |> 
    select(Cell_Type, cancer, Spearman_r) |> 
    pivot_wider(names_from = cancer, values_from = Spearman_r) |> 
    column_to_rownames("Cell_Type")
    
  pheatmap(
    as.matrix(plot_df),
    color = colorRampPalette(c("#2171b5", "white", "#de2d26"))(100),
    main = "Pan-Cancer: Immune-Risk Correlation Matrix",
    display_numbers = TRUE,
    fontsize_number = 8,
    number_color = "black",
    filename = file.path(pancancer_fig, "pancancer_immune_risk_cor_heatmap.png"),
    width = 8, height = 10
  )
  
  # ── Dotplot of Top Consistent Markers
  top_markers <- summary_cor |> head(15) |> pull(Cell_Type)
  dot_df <- all_cor |> filter(Cell_Type %in% top_markers)
  
  p_dot <- ggplot(dot_df, aes(x = cancer, y = Cell_Type, size = abs(Spearman_r), color = Spearman_r)) +
    geom_point() +
    scale_color_gradient2(low = "#2171b5", mid = "white", high = "#de2d26") +
    scale_size_continuous(range = c(2, 8)) +
    labs(
      title = "Top Consistent Immune Markers of Risk",
      subtitle = "Dot size = |Spearman r|; Color = Correlation Direction",
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  ggsave(file.path(pancancer_fig, "pancancer_immune_risk_dotplot.png"), 
         plot = p_dot, width = 8, height = 10, dpi = 600)

  message("  Pan-cancer immune aggregation complete.")
} else {
  message("  No per-cancer immune correlation files found.")
}
