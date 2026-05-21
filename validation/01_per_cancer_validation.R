source("validation/00_validation_functions.R")

for (cancer in CANCERS) {

  message("\n=== ", toupper(cancer), " ===")

  out_csv <- paste0("outputs/", cancer, "/validation/")
  out_fig <- paste0("outputs/", cancer, "/figures/validation/")
  dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

  # ── Load meta-analysis results
  meta <- tryCatch(load_meta(cancer), error = function(e) {
    message("  SKIPPING — cannot load meta: ", e$message); NULL
  })
  if (is.null(meta)) next

  # ── Download TCGA data via curatedTCGAData
  projects <- TCGA_CONFIG[[cancer]]
  message("  Projects: ", paste(projects, collapse = " + "))

  data <- if (length(projects) == 1) {
    mae <- download_tcga(projects)
    if (is.null(mae)) { message("  SKIPPING — download failed"); next }
    extract_mae_data(mae, projects)
  } else {
    combine_tcga(projects)
  }

  if (is.null(data)) { message("  SKIPPING — no data"); next }

  # Save
  dir.create(file.path("data/tcga_curated", cancer), showWarnings = FALSE, recursive = TRUE)
  saveRDS(data$expr, file.path("data/tcga_curated", cancer, paste0(cancer, "_expr.rds")))
  safe_export(data$meta, file.path("data/tcga_curated", cancer, paste0(cancer, "_meta.csv")))

  n_tumor  <- sum(data$meta$sample_type == "Tumor",  na.rm = TRUE)
  n_normal <- sum(data$meta$sample_type == "Normal", na.rm = TRUE)
  message("  Samples — Tumor: ", n_tumor, "  Normal: ", n_normal)

  # ── Limma differential expression
  message("  Running limma...")
  limma_res <- tryCatch(
    run_limma_validation(data$expr, data$meta, meta_df = meta),
    error = function(e) { message("  limma error: ", e$message); NULL }
  )
  if (is.null(limma_res)) next

  safe_export(limma_res$deseq, file.path(out_csv, paste0(cancer, "_TCGA_DEGs.csv")))

  # ── Compare with meta-analysis
  message("  Comparing with meta-analysis...")
  val <- tryCatch(
    compare_results(meta, limma_res$deseq),
    error = function(e) { message("  Comparison error: ", e$message); NULL }
  )

  if (is.null(val) || nrow(val) == 0) {
    message("  No overlapping genes for comparison")
    next
  }

  safe_export(val, file.path(out_csv, paste0(cancer, "_validation_results.csv")))

  # ── Metrics
  cor_res      <- cor.test(val$pooled_log2FC, val$log2FoldChange, method = "spearman")
  concord_rate <- round(mean(val$Concordant) * 100, 1)
  val_rate     <- round(mean(val$Validated)  * 100, 1)

  metrics <- tibble(
    cancer          = cancer,
    n_meta_sig      = sum(meta$rem_padj < PADJ_CUTOFF & abs(meta$pooled_log2FC) >= LFC_CUTOFF,
                          na.rm = TRUE),
    n_TCGA_DEGs     = nrow(limma_res$deseq |> filter(padj < PADJ_CUTOFF)),
    n_compared      = nrow(val),
    n_validated     = sum(val$Validated),
    concordance_pct = concord_rate,
    validation_pct  = val_rate,
    spearman_r      = round(cor_res$estimate, 4),
    spearman_p      = signif(cor_res$p.value, 4)
  )
  safe_export(metrics, file.path(out_csv, paste0(cancer, "_validation_metrics.csv")))

  message("  Compared: ", nrow(val), "  Validated: ", sum(val$Validated),
          "  Concordance: ", concord_rate, "%  r: ", round(cor_res$estimate, 3))

  # ── Figures
  save_correlation_plot(val, cancer,
    file.path(out_fig, paste0(cancer, "_correlation_scatter.png")))

  save_volcano_plot(limma_res$deseq, val, cancer,
    file.path(out_fig, paste0(cancer, "_volcano.png")))

  save_heatmap(limma_res$expr, limma_res$meta, val, cancer,
    file.path(out_fig, paste0(cancer, "_validated_heatmap.png")))

  message("  Done: ", cancer)
}

message("\nPer-cancer validation complete.")
