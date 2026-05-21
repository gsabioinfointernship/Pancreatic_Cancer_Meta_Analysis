source("survival_analysis/00_survival_functions.R")

for (cancer in CANCERS) {

  message("\n=== ", toupper(cancer), " ===")

  out_csv <- paste0("outputs/", cancer, "/survival/")
  out_fig <- paste0("outputs/", cancer, "/figures/survival/")
  dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

  # ── Load meta-analysis results
  meta <- tryCatch(load_meta(cancer), error = function(e) {
    message("  SKIPPING — cannot load meta: ", e$message); NULL
  })
  if (is.null(meta)) next

  # ── Load hub genes (from network analysis)
  hub_df <- load_hub_genes(cancer)
  if (is.null(hub_df)) {
    message("  Hub genes not found — using top sig DEGs")
  } else {
    message("  Hub genes loaded: ", nrow(hub_df))
  }

  # ── Select candidate genes
  candidates <- select_survival_genes(meta, hub_df, n = HUB_N)
  message("  Candidate genes for survival: ", nrow(candidates))

  if (nrow(candidates) == 0) {
    message("  SKIPPING — no candidate genes")
    next
  }

  # ── Download TCGA data via curatedTCGAData
  projects <- TCGA_CONFIG[[cancer]]
  message("  Projects: ", paste(projects, collapse = " + "))

  tcga_data <- if (length(projects) == 1) {
    mae <- download_tcga_survival(projects)
    if (is.null(mae)) { message("  SKIPPING — download failed"); next }
    extract_survival_data(mae, projects)
  } else {
    combine_tcga_survival(projects)
  }

  if (is.null(tcga_data)) { message("  SKIPPING — no TCGA data"); next }

  n_tumor  <- ncol(tcga_data$expr)
  n_clin   <- nrow(tcga_data$clinical)
  n_events <- sum(tcga_data$clinical$OS_status, na.rm = TRUE)
  message("  Tumor expression samples: ", n_tumor)
  message("  Patients with OS data: ", n_clin, "  Deaths: ", n_events)

  if (n_events < MIN_EVENTS) {
    message("  SKIPPING — too few death events (", n_events, " < ", MIN_EVENTS, ")")
    next
  }

  # Cache TCGA survival data
  cache_dir <- file.path("data/tcga_curated", cancer)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(tcga_data$expr,
          file.path(cache_dir, paste0(cancer, "_survival_expr.rds")))
  safe_export(tcga_data$clinical,
              file.path(cache_dir, paste0(cancer, "_survival_clinical.csv")))

  # ── Univariate Cox for all candidate genes
  message("  Running univariate Cox for ", nrow(candidates), " genes...")

  # Filter candidates to genes present in expression matrix
  query_genes <- intersect(candidates$Symbol, rownames(tcga_data$expr))
  message("  Genes in expression matrix: ", length(query_genes))

  if (length(query_genes) == 0) {
    message("  SKIPPING — no candidate genes in expression matrix")
    next
  }

  cox_res <- tryCatch(
    run_univariate_cox(tcga_data$expr, tcga_data$clinical, query_genes),
    error = function(e) { message("  Cox error: ", e$message); NULL }
  )

  if (is.null(cox_res) || nrow(cox_res) == 0) {
    message("  No Cox results (too few events per gene?)")
    next
  }

  # Annotate with meta-analysis direction
  cox_res <- cox_res |>
    left_join(meta |> select(Symbol, pooled_log2FC, rem_padj, n_studies),
              by = "Symbol") |>
    mutate(
      Direction    = ifelse(pooled_log2FC > 0, "Up", "Down"),
      Risk_aligned = case_when(
        Direction == "Up"   & HR > 1 ~ "Up-High risk",
        Direction == "Up"   & HR < 1 ~ "Up-Protective",
        Direction == "Down" & HR > 1 ~ "Down-High risk",
        Direction == "Down" & HR < 1 ~ "Down-Protective",
        TRUE ~ "Unknown"
      )
    )

  safe_export(cox_res,
              file.path(out_csv, paste0(cancer, "_univariate_cox.csv")))

  # Significant prognostic genes
  prog_sig <- cox_res |> filter(padj < COX_PADJ)
  message("  Prognostic genes (FDR < ", COX_PADJ, "): ", nrow(prog_sig))

  if (nrow(prog_sig) > 0) {
    safe_export(prog_sig,
                file.path(out_csv, paste0(cancer, "_prognostic_genes.csv")))
    message("  Top gene: ", prog_sig$Symbol[1],
            "  HR = ", round(prog_sig$HR[1], 3),
            "  p = ", signif(prog_sig$pvalue[1], 3))
  }

  # ── Summary metrics
  metrics <- tibble(
    cancer         = cancer,
    n_candidates   = nrow(candidates),
    n_tested       = nrow(cox_res),
    n_events       = n_events,
    n_prognostic   = nrow(prog_sig),
    n_high_risk    = sum(prog_sig$HR > 1, na.rm = TRUE),
    n_protective   = sum(prog_sig$HR <= 1, na.rm = TRUE),
    top_gene       = if (nrow(prog_sig) > 0) prog_sig$Symbol[1]          else NA_character_,
    top_HR         = if (nrow(prog_sig) > 0) round(prog_sig$HR[1], 3)    else NA_real_,
    top_HR_lo      = if (nrow(prog_sig) > 0) round(prog_sig$HR_lo[1], 3) else NA_real_,
    top_HR_hi      = if (nrow(prog_sig) > 0) round(prog_sig$HR_hi[1], 3) else NA_real_,
    top_pvalue     = if (nrow(prog_sig) > 0) signif(prog_sig$pvalue[1], 4) else NA_real_,
    top_padj       = if (nrow(prog_sig) > 0) signif(prog_sig$padj[1], 4)   else NA_real_
  )
  safe_export(metrics,
              file.path(out_csv, paste0(cancer, "_survival_metrics.csv")))

  # ── Figures
  message("  Saving figures...")

  # KM plots — top 9 significant prognostic genes
  if (nrow(prog_sig) > 0) {
    save_km_plots(
      expr_mat   = tcga_data$expr,
      clinical   = tcga_data$clinical,
      prog_genes = prog_sig$Symbol,
      cancer     = cancer,
      path       = file.path(out_fig, paste0(cancer, "_km_plots.png")),
      n_top      = 9
    )
    message("  KM plots saved")
  } else {
    message("  No significant genes for KM plots — saving all-gene forest plot only")
  }

  # Forest plot — all tested genes (significant highlighted)
  save_forest_plot(
    cox_df = cox_res,
    cancer = cancer,
    path   = file.path(out_fig, paste0(cancer, "_forest_plot.png"))
  )
  message("  Forest plot saved")

  message("  Done: ", cancer)
}

message("\nPer-cancer survival analysis complete.")
