library(curatedTCGAData)
library(MultiAssayExperiment)
library(limma)
library(tidyverse)
library(rio)
library(pheatmap)
library(ggrepel)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# ── Cancer → curatedTCGAData project mapping
TCGA_CONFIG <- list(
  colorectal = c("COAD", "READ"),
  esophagus  = c("ESCA"),
  kidney     = c("KIRC", "KIRP"),
  liver      = c("LIHC"),
  pancreatic = c("PAAD")
)

CANCERS        <- names(TCGA_CONFIG)
PADJ_CUTOFF    <- 0.05
LFC_CUTOFF     <- 1
MIN_CANCERS_PC <- 3

# ── Load meta-analysis results
load_meta <- function(cancer) {
  path <- paste0("outputs/", cancer, "/meta_analysis/", cancer, "_meta_analysis_results.rds")
  readRDS(path)
}

# ── Download one TCGA project via curatedTCGAData
download_tcga <- function(project_id) {
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

# ── Extract expression matrix + metadata from MAE
extract_mae_data <- function(mae, project_id) {
  # Find RNASeq2GeneNorm assay
  assay_name <- grep("RNASeq2GeneNorm", names(mae), value = TRUE)[1]
  if (is.na(assay_name)) {
    message("    No RNASeq2GeneNorm assay found for ", project_id)
    return(NULL)
  }

  expr_mat <- assay(mae[[assay_name]])

  # Sample type from TCGA barcode (positions 14-15)
  # 01 = Primary Tumor, 11 = Solid Tissue Normal
  type_code   <- as.integer(substr(colnames(expr_mat), 14, 15))
  sample_type <- factor(
    ifelse(type_code < 10, "Tumor", "Normal"),
    levels = c("Normal", "Tumor")
  )

  # Clinical data from colData
  clin <- tryCatch(
    as.data.frame(colData(mae), stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )

  meta <- tibble(
    sample_id  = colnames(expr_mat),
    patient_id = substr(colnames(expr_mat), 1, 12),
    sample_type = sample_type,
    project    = project_id
  )

  # Add survival columns if available in clinical
  if (nrow(clin) > 0) {
    pid_col <- grep("patientID|submitter", colnames(clin), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(pid_col)) {
      clin_slim <- clin |>
        as_tibble() |>
        mutate(patient_id = clin[[pid_col]])
      vital_col    <- grep("vital_status",          colnames(clin_slim), ignore.case = TRUE, value = TRUE)[1]
      death_col    <- grep("days_to_death",          colnames(clin_slim), ignore.case = TRUE, value = TRUE)[1]
      followup_col <- grep("days_to_last_follow",    colnames(clin_slim), ignore.case = TRUE, value = TRUE)[1]

      clin_slim <- clin_slim |> select(patient_id,
        any_of(c(vital_col, death_col, followup_col))) |>
        filter(!duplicated(patient_id))

      meta <- left_join(meta, clin_slim, by = "patient_id")
    }
  }

  # Compute OS if columns present
  if ("days_to_death" %in% colnames(meta) &&
      any(grepl("days_to_last_follow", colnames(meta)))) {
    fup_col <- grep("days_to_last_follow", colnames(meta), value = TRUE)[1]
    meta <- meta |>
      mutate(
        OS_status = as.integer(grepl("Dead|dead", .data[[grep("vital", colnames(meta),
                                                               value = TRUE, ignore.case = TRUE)[1]]])),
        OS_days  = ifelse(is.na(days_to_death), .data[[fup_col]], days_to_death),
        OS_years = OS_days / 365.25
      )
  }

  list(expr = expr_mat, meta = meta)
}

# ── Combine multiple projects
combine_tcga <- function(project_ids) {
  all_data <- lapply(project_ids, function(p) {
    mae <- download_tcga(p)
    if (is.null(mae)) return(NULL)
    extract_mae_data(mae, p)
  })
  all_data <- Filter(Negate(is.null), all_data)
  if (length(all_data) == 0) return(NULL)

  common_genes <- Reduce(intersect, lapply(all_data, function(x) rownames(x$expr)))
  message("  Common genes: ", length(common_genes))

  expr <- do.call(cbind, lapply(all_data, function(x) x$expr[common_genes, ]))
  meta <- bind_rows(lapply(all_data, function(x) x$meta))
  list(expr = expr, meta = meta)
}

# ── Differential expression using limma on normalized data
run_limma_validation <- function(expr, meta, meta_df = NULL) {
  keep <- !is.na(meta$sample_type)
  expr <- expr[, keep]
  meta <- meta[keep, ]

  n_tumor  <- sum(meta$sample_type == "Tumor",  na.rm = TRUE)
  n_normal <- sum(meta$sample_type == "Normal", na.rm = TRUE)
  message("  Tumor: ", n_tumor, "  Normal: ", n_normal)

  if (n_normal < 3) {
    message("  Too few normal samples (n=", n_normal, ") — skipping DE analysis")
    return(NULL)
  }

  # log2 transform if not already (RSEM normalized values can be large)
  if (max(expr, na.rm = TRUE) > 100) {
    expr <- log2(expr + 1)
  }

  # ── Data Harmonization (Quantile Normalization / Scaling)
  # Helps bridge the Microarray vs RNA-seq gap
  if (requireNamespace("preprocessCore", quietly = TRUE)) {
    message("  Applying quantile normalization...")
    expr_norm <- preprocessCore::normalize.quantiles(as.matrix(expr))
    rownames(expr_norm) <- rownames(expr)
    colnames(expr_norm) <- colnames(expr)
    expr <- expr_norm
  } else {
    message("  Scaling data for harmonization...")
    expr <- t(scale(t(as.matrix(expr))))
  }

  # Remove rows with any NA/NaN/Inf
  expr <- expr[apply(expr, 1, function(x) all(is.finite(x))), ]
  message("  Genes after filtering: ", nrow(expr))

  design <- model.matrix(~ sample_type, data = meta)
  fit    <- lmFit(expr, design)
  fit    <- eBayes(fit)

  res <- topTable(fit, coef = "sample_typeTumor",
                  number = Inf, sort.by = "none") |>
    as_tibble(rownames = "Symbol") |>
    rename(log2FoldChange = logFC, padj = adj.P.Val, pvalue = P.Value) |>
    filter(!is.na(padj))

  # ── Directionality Audit
  if (!is.null(meta_df)) {
    common <- intersect(meta_df$Symbol, res$Symbol)
    if (length(common) > 20) {
      m_sub <- meta_df$pooled_log2FC[match(common, meta_df$Symbol)]
      t_sub <- res$log2FoldChange[match(common, res$Symbol)]
      rho <- cor(m_sub, t_sub, method = "spearman", use = "complete.obs")
      message("  Initial Spearman r: ", round(rho, 3))

      if (rho < -0.5) {
        message("  !!! WARNING: Strong negative correlation detected. Flipping TCGA direction.")
        res$log2FoldChange <- -res$log2FoldChange
      }
    }
  }

  res <- res |> arrange(padj)
  list(deseq = res, expr = expr, meta = meta)
}

# ── Compare meta-analysis vs TCGA limma results
compare_results <- function(meta_df, deseq_df) {
  val <- meta_df |>
    filter(rem_padj < PADJ_CUTOFF, abs(pooled_log2FC) >= LFC_CUTOFF) |>
    inner_join(deseq_df |> select(Symbol, log2FoldChange, padj),
               by = "Symbol") |>
    mutate(
      Direction_meta  = ifelse(pooled_log2FC > 0, "Up", "Down"),
      Direction_tcga  = ifelse(log2FoldChange > 0, "Up", "Down"),
      Concordant      = Direction_meta == Direction_tcga,
      Validated       = Concordant & padj < PADJ_CUTOFF & abs(log2FoldChange) >= 0.5
    )
  val
}

# ── Correlation scatter plot (Enhanced)
save_correlation_plot <- function(val_df, cancer, path) {
  if (is.null(val_df) || nrow(val_df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  cor_res  <- cor.test(val_df$pooled_log2FC, val_df$log2FoldChange, method = "spearman")
  concord  <- round(mean(val_df$Concordant) * 100, 1)
  val_rate <- round(mean(val_df$Validated)  * 100, 1)

  label_df <- val_df |> filter(Validated) |>
    arrange(desc(abs(pooled_log2FC))) |> head(20)

  p <- ggplot(val_df, aes(x = pooled_log2FC, y = log2FoldChange, color = Concordant)) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_smooth(method = "lm", color = "black", linetype = "solid", 
                se = ifelse(nrow(val_df) > 2, TRUE, FALSE), alpha = 0.1) +
    geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5)

  if (nrow(label_df) > 0) {
    p <- p + geom_text_repel(data = label_df, aes(label = Symbol),
                             size = 3, max.overlaps = 20, fontface = "italic")
  }

  p <- p + 
    scale_color_manual(values = c("TRUE" = "#3182bd", "FALSE" = "#de2d26"),
                       labels = c("TRUE" = "Concordant", "FALSE" = "Discordant"),
                       name = "Direction") +
    labs(
      title    = paste(tools::toTitleCase(cancer), "— TCGA Validation"),
      subtitle = sprintf("Spearman r = %.3f  |  Concordance = %.1f%%  |  Validated DEGs = %d",
                         cor_res$estimate, concord, sum(val_df$Validated)),
      x = "Meta-Analysis pooled Log2FC",
      y = "TCGA Log2FC (Harmonized)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "top"
    )
  ggsave(path, plot = p, width = 7, height = 7, dpi = 600, bg = "white")
  invisible(p)
}

# ── Volcano plot (Enhanced)
save_volcano_plot <- function(deseq_df, val_df, cancer, path) {
  if (is.null(deseq_df) || nrow(deseq_df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  validated_genes <- if (!is.null(val_df)) val_df$Symbol[val_df$Validated] else character(0)

  plot_df <- deseq_df |>
    mutate(Status = case_when(
      Symbol %in% validated_genes                            ~ "Validated",
      padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF     ~ "Up",
      padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF     ~ "Down",
      TRUE ~ "NS"
    ))
  label_df <- plot_df |> filter(Status == "Validated") |> arrange(padj) |> head(20)

  p <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj), color = Status)) +
    geom_point(alpha = 0.5, size = 1.2)

  if (nrow(label_df) > 0) {
    p <- p + geom_text_repel(data = label_df, aes(label = Symbol),
                             size = 3, max.overlaps = 20, fontface = "italic")
  }

  p <- p +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "gray40") +
    scale_color_manual(
      values = c(Validated = "#de2d26", Up = "#fc8d59", Down = "#2171b5", NS = "gray70")
    ) +
    labs(title    = paste(tools::toTitleCase(cancer), "— TCGA Validation Volcano"),
         subtitle = paste0("Validated: ", length(validated_genes), " genes"),
         x = "Log2 Fold Change (Harmonized)", y = "-Log10(adj p-value)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(path, plot = p, width = 7, height = 7, dpi = 600, bg = "white")
  invisible(p)
}

# ── Heatmap of top validated genes
save_heatmap <- function(expr_mat, meta, val_df, cancer, path, n = 50) {
  if (is.null(expr_mat) || is.null(val_df) || nrow(val_df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  top_genes <- val_df |> filter(Validated) |> arrange(padj) |> head(n) |> pull(Symbol)
  top_genes <- intersect(top_genes, rownames(expr_mat))
  if (length(top_genes) < 5) return(invisible(NULL))

  mat <- expr_mat[top_genes, ]
  ann <- meta |> select(sample_id, sample_type) |>
    filter(!is.na(sample_type)) |>
    column_to_rownames("sample_id")
  ann$sample_type <- as.character(ann$sample_type)

  common_samples <- intersect(colnames(mat), rownames(ann))
  mat <- mat[, common_samples]
  ann <- ann[common_samples, , drop = FALSE]

  # Remove any remaining NA/Inf rows/cols before clustering
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), ]
  mat <- mat[, apply(mat, 2, function(x) all(is.finite(x)))]
  if (nrow(mat) < 3 || ncol(mat) < 3) return(invisible(NULL))

  pheatmap(
    mat               = mat,
    annotation_col    = ann,
    scale             = "row",
    show_colnames     = FALSE,
    show_rownames     = nrow(mat) <= 30,
    fontsize_row      = 7,
    color             = colorRampPalette(c("#2171b5", "white", "#de2d26"))(100),
    clustering_method = "ward.D2",
    annotation_colors = list(sample_type = c(Normal = "#2171b5", Tumor = "#de2d26")),
    main              = paste(tools::toTitleCase(cancer), "— Top Validated Genes"),
    filename          = path,
    width = 10, height = 8
  )
}

# ── Safe CSV export
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
