source("machine_learning/00_ml_functions.R")

for (cancer in CANCERS) {

  message("\n=== ", toupper(cancer), " ===")

  # 1. Setup directories
  out_csv <- paste0("outputs/", cancer, "/machine_learning/")
  out_fig <- paste0("outputs/", cancer, "/figures/machine_learning/")
  dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

  # 2. Load meta and hub genes
  meta <- tryCatch(load_meta(cancer), error = function(e) {
    message("  SKIPPING — cannot load meta: ", e$message); NULL
  })
  if (is.null(meta)) next

  hub_df <- load_hub_genes(cancer)
  if (is.null(hub_df)) {
    message("  Hub genes not found — using top sig DEGs as features")
  }

  # 3. Select candidate genes for ML
  candidates <- select_ml_genes(meta, hub_df, n = HUB_N)
  message("  ML candidate genes: ", nrow(candidates))
  if (nrow(candidates) == 0) {
    message("  SKIPPING — no candidate genes")
    next
  }

  # 4. Download TCGA data
  projects <- TCGA_CONFIG[[cancer]]
  message("  TCGA Projects: ", paste(projects, collapse = " + "))

  tcga_data <- combine_tcga_ml(projects)
  if (is.null(tcga_data)) {
    message("  SKIPPING — no TCGA data")
    next
  }

  # ============================================================================
  # LASSO-COX PROGNOSTIC MODEL
  # ============================================================================
  message("\n  --- LASSO-COX PROGNOSTIC MODEL ---")

  lasso_data <- prepare_lasso_data(
    expr_mat    = tcga_data$expr,
    sample_type = tcga_data$sample_type,
    clinical    = tcga_data$clinical,
    genes       = candidates$Symbol
  )

  lasso_res <- NULL
  risk_eval <- NULL

  if (!is.null(lasso_data)) {
    lasso_res <- run_lasso_cox(lasso_data)

    if (!is.null(lasso_res)) {
      risk_eval <- compute_and_evaluate_risk(lasso_res, lasso_data, cancer)

      # Export results
      coef_df <- as.data.frame(as.matrix(lasso_res$coefs)) |>
        rownames_to_column("Symbol") |>
        rename_with(~ "Coefficient", 2)

      safe_export(coef_df,     file.path(out_csv, paste0(cancer, "_lasso_coefficients.csv")))
      safe_export(risk_eval$df, file.path(out_csv, paste0(cancer, "_patient_risk_scores.csv")))

      # Plots
      save_lasso_cv_plot(lasso_res$cv_fit, cancer,
                         file.path(out_fig, paste0(cancer, "_lasso_cv_curve.png")))
      save_risk_km_plot(risk_eval$df, cancer, risk_eval$c_index, risk_eval$km_pval,
                        file.path(out_fig, paste0(cancer, "_risk_km_plot.png")))
      save_risk_distribution(risk_eval$df, cancer,
                             file.path(out_fig, paste0(cancer, "_risk_distribution.png")))
      save_coef_plot(coef_df, cancer,
                     file.path(out_fig, paste0(cancer, "_lasso_coef_barplot.png")))

      message("  LASSO-Cox complete.")
    } else {
      message("  LASSO-Cox failed to converge or select genes.")
    }
  }

  # ============================================================================
  # RANDOM FOREST DIAGNOSTIC MODEL
  # ============================================================================
  message("\n  --- RANDOM FOREST DIAGNOSTIC MODEL ---")

  rf_data <- prepare_rf_data(
    expr_mat    = tcga_data$expr,
    sample_type = tcga_data$sample_type,
    genes       = candidates$Symbol
  )

  rf_res <- NULL

  if (!is.null(rf_data)) {
    rf_res <- run_rf_classifier(rf_data)

    if (!is.null(rf_res)) {
      # Export results
      safe_export(rf_res$importance,
                  file.path(out_csv, paste0(cancer, "_rf_importance.csv")))

      # Metrics table
      rf_metrics <- data.frame(
        Metric = names(rf_res$cm$overall),
        Value  = as.numeric(rf_res$cm$overall)
      ) |>
        bind_rows(data.frame(
          Metric = names(rf_res$cm$byClass),
          Value  = as.numeric(rf_res$cm$byClass)
        )) |>
        bind_rows(data.frame(Metric = "AUC", Value = rf_res$auc))

      safe_export(rf_metrics, file.path(out_csv, paste0(cancer, "_rf_metrics.csv")))

      # Plots
      save_roc_plot(rf_res, cancer,
                    file.path(out_fig, paste0(cancer, "_rf_roc_curve.png")))
      save_importance_plot(rf_res, cancer,
                           file.path(out_fig, paste0(cancer, "_rf_importance_barplot.png")))

      message("  Random Forest complete.")
    } else {
      message("  Random Forest failed.")
    }
  }

  # ============================================================================
  # SUMMARY METRICS
  # ============================================================================
  summary_row <- tibble(
    cancer       = cancer,
    lasso_genes  = if (!is.null(lasso_res)) lasso_res$n_selected else NA_integer_,
    lasso_cindex = if (!is.null(risk_eval)) round(risk_eval$c_index, 3) else NA_real_,
    lasso_km_p   = if (!is.null(risk_eval)) signif(risk_eval$km_pval, 3) else NA_real_,
    rf_auc       = if (!is.null(rf_res)) round(rf_res$auc, 3) else NA_real_,
    rf_accuracy  = if (!is.null(rf_res)) round(rf_res$accuracy, 3) else NA_real_,
    rf_sensitivity = if (!is.null(rf_res)) round(rf_res$sens, 3) else NA_real_,
    rf_specificity = if (!is.null(rf_res)) round(rf_res$spec, 3) else NA_real_
  )

  safe_export(summary_row, file.path(out_csv, paste0(cancer, "_ml_summary.csv")))

  message("\nDone: ", cancer)
}

message("\nPer-cancer machine learning complete.")
