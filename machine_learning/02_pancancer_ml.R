source("machine_learning/00_ml_functions.R")

# 1. Setup pan-cancer output directory
pancancer_csv <- "outputs/pancancer/machine_learning/"
pancancer_fig <- "outputs/pancancer/figures/machine_learning/"
dir.create(pancancer_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(pancancer_fig, showWarnings = FALSE, recursive = TRUE)

message("\n=== PAN-CANCER MACHINE LEARNING AGGREGATION ===")

# ============================================================================
# 1. AGGREGATE SUMMARY METRICS
# ============================================================================
message("  Aggregating summary metrics...")
summary_files <- list.files("outputs", pattern = "_ml_summary.csv",
                            recursive = TRUE, full.names = TRUE)

if (length(summary_files) > 0) {
  all_metrics <- lapply(summary_files, read.csv) |>
    bind_rows() |>
    arrange(desc(rf_auc))

  safe_export(all_metrics, file.path(pancancer_csv, "pancancer_ml_summary_metrics.csv"))
  message("  Saved: pancancer_ml_summary_metrics.csv")
} else {
  message("  No summary metrics found.")
}

# ============================================================================
# 2. AGGREGATE LASSO-COX COEFFICIENTS
# ============================================================================
message("  Aggregating LASSO-Cox coefficients...")
coef_files <- list.files("outputs", pattern = "_lasso_coefficients.csv",
                         recursive = TRUE, full.names = TRUE)

if (length(coef_files) > 0) {
  all_coefs <- lapply(coef_files, function(f) {
    cancer_name <- gsub("_lasso_coefficients.csv", "", basename(f))
    read.csv(f) |> mutate(cancer = cancer_name)
  }) |> bind_rows()

  # Identify consistently selected genes
  consistent_prognostic <- all_coefs |>
    group_by(Symbol) |>
    summarise(
      n_cancers    = n_distinct(cancer),
      cancers      = paste(unique(cancer), collapse = "; "),
      mean_coef    = mean(Coefficient),
      median_coef  = median(Coefficient),
      .groups = "drop"
    ) |>
    arrange(desc(n_cancers), desc(abs(mean_coef)))

  safe_export(all_coefs, file.path(pancancer_csv, "pancancer_all_lasso_coefficients.csv"))
  safe_export(consistent_prognostic,
              file.path(pancancer_csv, "pancancer_consistent_prognostic_features.csv"))

  # Top consistent prognostic features table
  message("  Consistent prognostic genes (selected in >= 2 cancers): ",
          sum(consistent_prognostic$n_cancers >= 2))
} else {
  message("  No LASSO coefficients found.")
}

# ============================================================================
# 3. AGGREGATE RF FEATURE IMPORTANCE
# ============================================================================
message("  Aggregating Random Forest feature importance...")
rf_imp_files <- list.files("outputs", pattern = "_rf_importance.csv",
                           recursive = TRUE, full.names = TRUE)

if (length(rf_imp_files) > 0) {
  all_rf_imp <- lapply(rf_imp_files, function(f) {
    cancer_name <- gsub("_rf_importance.csv", "", basename(f))
    read.csv(f) |> mutate(cancer = cancer_name)
  }) |> bind_rows()

  # Rank features by mean importance across cancers
  top_diagnostic <- all_rf_imp |>
    group_by(Symbol) |>
    summarise(
      n_cancers       = n_distinct(cancer),
      cancers         = paste(unique(cancer), collapse = "; "),
      mean_gini_dec   = mean(MeanDecreaseGini, na.rm = TRUE),
      median_gini_dec = median(MeanDecreaseGini, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(n_cancers), desc(mean_gini_dec))

  safe_export(all_rf_imp, file.path(pancancer_csv, "pancancer_all_rf_importance.csv"))
  safe_export(top_diagnostic,
              file.path(pancancer_csv, "pancancer_top_diagnostic_features.csv"))

  message("  Top diagnostic features saved.")
} else {
  message("  No RF importance found.")
}

# ============================================================================
# 4. PAN-CANCER VISUALIZATIONS
# ============================================================================
message("  Generating pan-cancer ML plots...")

# A. Diagnostic Performance Dotplot
if (exists("all_metrics")) {
  p_auc <- ggplot(all_metrics, aes(x = reorder(cancer, rf_auc), y = rf_auc, fill = rf_auc)) +
    geom_bar(stat = "identity", show.legend = FALSE, width = 0.7) +
    geom_text(aes(label = round(rf_auc, 3)), vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_gradient(low = "#fee0d2", high = "#de2d26") +
    labs(title = "Diagnostic Performance Across Cancers",
         subtitle = "Random Forest Classifier AUC (Hold-out Test Set)",
         x = NULL, y = "Area Under the ROC Curve (AUC)") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(face = "bold")
    )

  ggsave(file.path(pancancer_fig, "pancancer_rf_auc_barplot.png"),
         plot = p_auc, width = 7, height = 5, dpi = 600, bg = "white")

  p_cindex <- ggplot(all_metrics, aes(x = reorder(cancer, lasso_cindex), y = lasso_cindex, fill = lasso_cindex)) +
    geom_bar(stat = "identity", show.legend = FALSE, width = 0.7) +
    geom_text(aes(label = round(lasso_cindex, 3)), vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
    labs(title = "Prognostic Performance Across Cancers",
         subtitle = "LASSO-Cox Risk Model Concordance Index",
         x = NULL, y = "C-index") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(face = "bold")
    )

  ggsave(file.path(pancancer_fig, "pancancer_lasso_cindex_barplot.png"),
         plot = p_cindex, width = 7, height = 5, dpi = 600, bg = "white")
}

# B. Top Consistent Prognostic Features Heatmap
if (exists("consistent_prognostic") && nrow(consistent_prognostic) > 0) {
  top_prog_genes <- consistent_prognostic |>
    filter(n_cancers >= 2) |>
    head(30) |>
    pull(Symbol)

  if (length(top_prog_genes) > 2) {
    plot_df <- all_coefs |>
      filter(Symbol %in% top_prog_genes) |>
      mutate(Symbol = factor(Symbol, levels = rev(top_prog_genes)))

    p_heat <- ggplot(plot_df, aes(x = cancer, y = Symbol, fill = Coefficient)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient2(low = "#3182bd", mid = "white", high = "#de2d26",
                           midpoint = 0, name = "LASSO\nCoef") +
      labs(title = "Consistent Prognostic Features",
           subtitle = "Features selected in ≥ 2 cancers",
           x = NULL, y = NULL) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "italic"),
        panel.grid = element_blank()
      )

    ggsave(file.path(pancancer_fig, "pancancer_consistent_prognostic_heatmap.png"),
           plot = p_heat, width = 6, height = 8, dpi = 600, bg = "white")
  }
}

message("\nPan-cancer ML aggregation complete.")
