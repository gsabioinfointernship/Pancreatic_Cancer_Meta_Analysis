library(survival)
library(survminer)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(gridExtra)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# ── Cancer → TCGA project mapping (same as validation)
TCGA_CONFIG <- list(
  colorectal = c("COAD", "READ"),
  esophagus  = c("ESCA"),
  kidney     = c("KIRC", "KIRP"),
  liver      = c("LIHC"),
  pancreatic = c("PAAD")
)

CANCERS      <- names(TCGA_CONFIG)
PADJ_CUTOFF  <- 0.05
LFC_CUTOFF   <- 1
HUB_N        <- 100    # top N hub genes to test per cancer
COX_PADJ     <- 0.05   # FDR cutoff for prognostic significance
MIN_EVENTS   <- 10     # minimum deaths required for reliable KM/Cox
MIN_CANCERS  <- 2      # pan-cancer: gene must be prognostic in >= N cancers

# ── Load meta-analysis results
load_meta <- function(cancer) {
  path <- paste0("outputs/", cancer, "/meta_analysis/", cancer, "_meta_analysis_results.rds")
  readRDS(path)
}

# ── Load hub genes from network analysis
load_hub_genes <- function(cancer) {
  path <- paste0("outputs/", cancer, "/network/", cancer, "_hub_genes.csv")
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}

# ── Download TCGA project via curatedTCGAData
download_tcga_survival <- function(project_id) {
  message("    Downloading: ", project_id)
  mae <- tryCatch(
    curatedTCGAData(
      diseaseCode = project_id,
      assays      = "RNASeq2GeneNorm",
      version     = "2.1.1",
      dry.run     = FALSE
    ),
    error = function(e) { message("    Error: ", e$message); NULL }
  )
  mae
}

# ── Extract tumor-only expression + OS clinical from MAE
extract_survival_data <- function(mae, project_id) {
  assay_name <- grep("RNASeq2GeneNorm", names(mae), value = TRUE)[1]
  if (is.na(assay_name)) {
    message("    No RNASeq2GeneNorm assay for ", project_id)
    return(NULL)
  }

  expr_mat  <- assay(mae[[assay_name]])

  # Keep PRIMARY TUMOR samples only (barcode position 14-15 < 10)
  type_code <- as.integer(substr(colnames(expr_mat), 14, 15))
  expr_mat  <- expr_mat[, type_code < 10, drop = FALSE]

  if (ncol(expr_mat) == 0) {
    message("    No tumor samples for ", project_id)
    return(NULL)
  }

  # Clinical data
  clin <- tryCatch(
    as.data.frame(colData(mae), stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )

  if (nrow(clin) == 0) {
    message("    No clinical data for ", project_id)
    return(NULL)
  }

  # Identify key columns
  pid_col      <- grep("patientID|submitter_id", colnames(clin), ignore.case = TRUE, value = TRUE)[1]
  vital_col    <- grep("vital_status",           colnames(clin), ignore.case = TRUE, value = TRUE)[1]
  death_col    <- grep("days_to_death",           colnames(clin), ignore.case = TRUE, value = TRUE)[1]
  followup_col <- grep("days_to_last_follow",     colnames(clin), ignore.case = TRUE, value = TRUE)[1]

  if (is.na(pid_col) || is.na(vital_col)) {
    message("    Missing clinical columns (pid or vital_status) for ", project_id)
    return(NULL)
  }

  clin_slim <- clin |>
    as_tibble() |>
    mutate(patient_id = .data[[pid_col]]) |>
    select(patient_id, any_of(c(vital_col, death_col, followup_col))) |>
    filter(!duplicated(patient_id))

  # Build OS
  clin_slim <- clin_slim |>
    mutate(
      OS_status = as.integer(tolower(.data[[vital_col]]) %in% c("dead", "deceased", "1")),
      OS_days   = {
        d <- if (!is.na(death_col) && death_col %in% colnames(clin_slim))
               as.numeric(.data[[death_col]]) else NA_real_
        f <- if (!is.na(followup_col) && followup_col %in% colnames(clin_slim))
               as.numeric(.data[[followup_col]]) else NA_real_
        dplyr::coalesce(d, f)
      }
    ) |>
    filter(!is.na(OS_days), is.finite(OS_days), OS_days > 0) |>
    mutate(OS_years = OS_days / 365.25)

  n_events <- sum(clin_slim$OS_status, na.rm = TRUE)
  message("    ", project_id, ": ", ncol(expr_mat), " tumor samples, ",
          nrow(clin_slim), " patients with OS, ", n_events, " deaths")

  list(expr = expr_mat, clinical = clin_slim)
}

# ── Combine multiple TCGA projects
combine_tcga_survival <- function(project_ids) {
  all_data <- lapply(project_ids, function(p) {
    mae <- download_tcga_survival(p)
    if (is.null(mae)) return(NULL)
    extract_survival_data(mae, p)
  })
  all_data <- Filter(Negate(is.null), all_data)
  if (length(all_data) == 0) return(NULL)

  common_genes <- Reduce(intersect, lapply(all_data, function(x) rownames(x$expr)))
  message("  Common genes across projects: ", length(common_genes))

  expr <- do.call(cbind, lapply(all_data, function(x) x$expr[common_genes, ]))
  clin <- bind_rows(lapply(all_data, function(x) x$clinical)) |>
    filter(!duplicated(patient_id))

  list(expr = expr, clinical = clin)
}

# ── Select candidate genes: top N hub genes that are significant in meta-analysis
select_survival_genes <- function(meta_df, hub_df, n = HUB_N) {
  sig_meta <- meta_df |>
    filter(rem_padj < PADJ_CUTOFF, abs(pooled_log2FC) >= LFC_CUTOFF) |>
    select(Symbol, pooled_log2FC, rem_padj)

  if (!is.null(hub_df) && nrow(hub_df) > 0) {
    candidates <- hub_df |>
      inner_join(sig_meta, by = "Symbol") |>
      arrange(desc(Degree)) |>
      head(n)
  } else {
    # Fallback: top sig DEGs by adjusted p-value
    candidates <- sig_meta |>
      arrange(rem_padj) |>
      head(n)
  }
  candidates
}

# ── Build per-gene survival dataframe (expression + OS, one row per patient)
build_gene_surv_df <- function(expr_mat, clinical, gene) {
  if (!gene %in% rownames(expr_mat)) return(NULL)

  gene_expr <- as.numeric(expr_mat[gene, ])
  names(gene_expr) <- colnames(expr_mat)

  # log2 transform if not already done
  if (max(gene_expr, na.rm = TRUE) > 100) gene_expr <- log2(gene_expr + 1)

  df <- tibble(
    sample_id  = names(gene_expr),
    patient_id = substr(names(gene_expr), 1, 12),
    expression = gene_expr
  ) |>
    filter(!is.na(expression), is.finite(expression)) |>
    left_join(clinical |> select(patient_id, OS_days, OS_status, OS_years),
              by = "patient_id") |>
    filter(!is.na(OS_days), !is.na(OS_status), OS_days > 0) |>
    filter(!duplicated(patient_id))  # one sample per patient

  if (nrow(df) < 10) return(NULL)

  # Median split
  med        <- median(df$expression, na.rm = TRUE)
  df$group   <- factor(ifelse(df$expression >= med, "High", "Low"), levels = c("Low", "High"))
  df
}

# ── Univariate Cox regression for a list of genes (continuous expression)
run_univariate_cox <- function(expr_mat, clinical, genes) {
  results <- lapply(genes, function(gene) {
    df <- build_gene_surv_df(expr_mat, clinical, gene)
    if (is.null(df) || nrow(df) < 10) return(NULL)
    if (sum(df$OS_status, na.rm = TRUE) < MIN_EVENTS) return(NULL)

    cx <- tryCatch(
      coxph(Surv(OS_days, OS_status) ~ expression, data = df),
      error = function(e) NULL
    )
    if (is.null(cx)) return(NULL)

    s <- summary(cx)
    tibble(
      Symbol   = gene,
      n        = nrow(df),
      n_events = sum(df$OS_status, na.rm = TRUE),
      HR       = s$conf.int[1, "exp(coef)"],
      HR_lo    = s$conf.int[1, "lower .95"],
      HR_hi    = s$conf.int[1, "upper .95"],
      pvalue   = s$coefficients[1, "Pr(>|z|)"],
      z        = s$coefficients[1, "z"]
    )
  })

  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(NULL)

  bind_rows(results) |>
    mutate(padj = p.adjust(pvalue, method = "BH")) |>
    arrange(pvalue)
}

# ── KM plots grid for top prognostic genes
save_km_plots <- function(expr_mat, clinical, prog_genes, cancer, path, n_top = 9) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  genes_to_plot <- head(prog_genes, n_top)
  if (length(genes_to_plot) == 0) return(invisible(NULL))

  plot_list <- lapply(genes_to_plot, function(gene) {
    df <- build_gene_surv_df(expr_mat, clinical, gene)
    if (is.null(df)) return(NULL)

    fit <- tryCatch(
      survfit(Surv(OS_days, OS_status) ~ group, data = df),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    lrt   <- tryCatch(survdiff(Surv(OS_days, OS_status) ~ group, data = df), error = function(e) NULL)
    p_val <- if (!is.null(lrt)) 1 - pchisq(lrt$chisq, df = length(lrt$n) - 1) else NA

    ggsurvplot(
      fit, data = df,
      palette     = c("#2171b5", "#de2d26"),
      title       = gene,
      subtitle    = if (!is.na(p_val)) sprintf("Log-rank p = %.3g", p_val) else "",
      xlab        = "Days",
      legend.labs = c("Low", "High"),
      ggtheme     = theme_bw(base_size = 9),
      risk.table  = FALSE,
      conf.int    = FALSE
    )$plot
  })

  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 0) return(invisible(NULL))

  nc  <- min(3, length(plot_list))
  nr  <- ceiling(length(plot_list) / nc)
  ht  <- max(4, nr * 4)
  wd  <- nc * 4

  tryCatch({
    png(path, width = wd, height = ht, units = "in", res = 150)
    do.call(grid.arrange, c(plot_list, list(ncol = nc)))
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("  KM grid save error: ", e$message)
  })
  invisible(NULL)
}

# ── Forest plot of univariate HR (top significant genes)
save_forest_plot <- function(cox_df, cancer, path, n_top = 30) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (is.null(cox_df) || nrow(cox_df) == 0) return(invisible(NULL))

  plot_df <- cox_df |>
    filter(padj < COX_PADJ) |>
    arrange(pvalue) |>
    head(n_top) |>
    mutate(
      Symbol = factor(Symbol, levels = rev(Symbol)),
      Risk   = ifelse(HR > 1, "High risk", "Protective")
    )

  if (nrow(plot_df) == 0) {
    # Show top 15 even if not significant
    plot_df <- cox_df |>
      arrange(pvalue) |>
      head(15) |>
      mutate(
        Symbol = factor(Symbol, levels = rev(Symbol)),
        Risk   = ifelse(HR > 1, "High risk", "Protective")
      )
  }

  p <- ggplot(plot_df, aes(x = HR, y = Symbol, color = Risk)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = HR_lo, xmax = HR_hi), height = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c("High risk" = "#de2d26", "Protective" = "#2171b5")) +
    scale_x_log10() +
    labs(
      title = paste(tools::toTitleCase(cancer), "— Prognostic Genes (Univariate Cox)"),
      x     = "Hazard Ratio (log scale, 95% CI)",
      y     = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )

  ht <- max(5, nrow(plot_df) * 0.35 + 2)
  ggsave(path, plot = p, width = 8, height = min(20, ht), dpi = 300)
  invisible(p)
}

# ── Pan-cancer HR heatmap (log2 HR, rows = genes, cols = cancers)
save_hr_heatmap <- function(hr_matrix, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (is.null(hr_matrix) || nrow(hr_matrix) < 3 || ncol(hr_matrix) < 2) return(invisible(NULL))

  mat     <- as.matrix(hr_matrix)
  mat_log <- log2(mat)  # HR=1 → 0, HR>1 → positive, HR<1 → negative

  # Remove rows that are all NA
  mat_log <- mat_log[rowSums(!is.na(mat_log)) >= 2, , drop = FALSE]
  if (nrow(mat_log) < 3) return(invisible(NULL))

  ht <- max(6, nrow(mat_log) * 0.25 + 3)

  pheatmap(
    mat_log,
    color             = colorRampPalette(c("#2171b5", "white", "#de2d26"))(100),
    na_col            = "gray90",
    cluster_rows      = TRUE,
    cluster_cols      = nrow(mat_log) > 3,
    clustering_method = "ward.D2",
    show_rownames     = nrow(mat_log) <= 50,
    fontsize_row      = 8,
    fontsize_col      = 10,
    main              = "Prognostic Hub Genes — Log2(Hazard Ratio) Across Cancers",
    filename          = path,
    width  = 8,
    height = min(20, ht)
  )
}

# ── Safe CSV export (same pattern as other modules)
safe_export <- function(df, path) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- sapply(df[[nm]], function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        tryCatch(paste(unique(rapply(x, as.character, how = "unlist")), collapse = "; "),
                 error = function(e) NA_character_)
      }, USE.NAMES = FALSE)
    }
  }
  write.csv(df, path, row.names = FALSE)
}
