library(GSVA)
library(tidyverse)
library(rio)
library(ggplot2)
library(pheatmap)
library(corrplot)
library(ggpubr)

# ── Shared Configuration
CANCERS <- c("colorectal", "esophagus", "kidney", "liver", "pancreatic")

# ── Immune Cell Marker Signatures (Bindea et al., 2013 / Charoentong et al., 2017)
# For high-impact, these 28 cell types are the gold standard.
IMMUNE_SIGNATURES <- list(
  Activated_B_cell = c("CD79A", "CD79B", "CD19", "MS4A1", "CD22", "CD24", "CD72", "BLNK", "PAX5", "SYK"),
  Activated_CD4_T_cell = c("CD4", "IL2RA", "CD25", "CD69", "CD44", "TNFRSF4", "TNFRSF18", "ICOS", "TNFRSF9"),
  Activated_CD8_T_cell = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "IFNG", "CD69", "TNFRSF9"),
  Activated_dendritic_cell = c("CD80", "CD83", "CD86", "LAMP3", "CCL19", "CCL21", "CCR7"),
  B_cell = c("CD19", "MS4A1", "CD79A", "CD79B", "CD22", "CD24", "CD72", "CD138"),
  CD56bright_natural_killer_cell = c("NCAM1", "NCR3", "KLRC1", "KLRD1", "KLRB1"),
  CD56dim_natural_killer_cell = c("NCAM1", "FCGR3A", "FGFBP2", "CX3CR1"),
  Central_memory_CD4_T_cell = c("CD4", "CCR7", "SELL", "CD27", "CD28", "IL7R"),
  Central_memory_CD8_T_cell = c("CD8A", "CD8B", "CCR7", "SELL", "CD27", "CD28", "IL7R"),
  Effector_memory_CD4_T_cell = c("CD4", "CCR5", "CXCR3", "SELL", "CD44", "FASL"),
  Effector_memory_CD8_T_cell = c("CD8A", "CD8B", "CCR5", "CXCR3", "SELL", "CD44", "FASL"),
  Eosinophil = c("CCR3", "SIGLEC8", "IL5RA", "PRG2", "RNASE3", "EPX"),
  Gamma_delta_T_cell = c("TRGC1", "TRGC2", "TRDC", "TRGV9", "TRDV2"),
  Immature_B_cell = c("CD19", "CD24", "CD38", "MS4A1"),
  Immature_dendritic_cell = c("CD1C", "CD1A", "CD209", "CD141"),
  MDSC = c("CD33", "ITGAM", "CD14", "FUT4", "S100A8", "S100A9", "ARG1", "NOS2"),
  Macrophage = c("CD68", "CD163", "CD14", "CSF1R", "MAFB"),
  Mast_cell = c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT", "FCER1A"),
  Memory_B_cell = c("CD19", "CD27", "CD38", "MS4A1", "IGHG1", "IGHA1"),
  Monocyte = c("CD14", "FCGR3A", "ITGAM", "CCR2", "CD33"),
  Natural_killer_cell = c("NCAM1", "NCR1", "NCR3", "KLRB1", "KLRD1", "KLRF1"),
  Natural_killer_T_cell = c("CD3D", "CD3E", "CD3G", "NCAM1", "NCR1", "KLRB1"),
  Neutrophil = c("FCGR3B", "CXCR1", "CXCR2", "S100A8", "S100A9", "MMP9"),
  Plasmacytoid_dendritic_cell = c("CLEC4C", "IL3RA", "NRP1", "TCF4", "LILRA4"),
  Regulatory_T_cell = c("CD4", "FOXP3", "IL2RA", "CD25", "IKZF2"),
  T_follicular_helper_cell = c("CD4", "CXCR5", "BCL6", "PDCD1", "ICOS", "IL21"),
  Type_1_T_helper_cell = c("CD4", "TBX21", "IFNG", "IL12RB2", "STAT1"),
  Type_17_T_helper_cell = c("CD4", "RORC", "IL17A", "IL17F", "IL22", "CCR6"),
  Type_2_T_helper_cell = c("CD4", "GATA3", "IL4", "IL5", "IL13", "STAT6")
)

# ── Helper: Load TCGA Expression
load_tcga_expr <- function(cancer) {
  path <- file.path("data/tcga_curated", cancer, paste0(cancer, "_expr.rds"))
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

# ── Helper: Load ML Risk Scores
load_risk_scores <- function(cancer) {
  path <- file.path("outputs", cancer, "machine_learning", paste0(cancer, "_patient_risk_scores.csv"))
  if (!file.exists(path)) return(NULL)
  read.csv(path)
}

# ── Core: Calculate ssGSEA Scores
run_immune_ssgsea <- function(expr_mat) {
  message("    Running ssGSEA for 29 immune cell types...")
  
  # Ensure matrix format
  if (!is.matrix(expr_mat)) expr_mat <- as.matrix(expr_mat)
  
  # Remove rows with any NAs (critical for clean GSVA results)
  na_rows <- rowSums(is.na(expr_mat)) > 0
  if (any(na_rows)) {
    message("    Removing ", sum(na_rows), " genes with NA values...")
    expr_mat <- expr_mat[!na_rows, ]
  }
  
  # log2 transform if needed
  if (max(expr_mat, na.rm = TRUE) > 100) {
    message("    Log-transforming expression data...")
    expr_mat <- log2(expr_mat + 1)
  }
  
  # Filter signatures to genes present in data
  present_genes <- rownames(expr_mat)
  filtered_sigs <- lapply(IMMUNE_SIGNATURES, function(x) intersect(x, present_genes))
  
  # Check if any signatures have too few genes
  sig_counts <- sapply(filtered_sigs, length)
  if (all(sig_counts == 0)) {
    message("    ERROR: No signature genes found in expression matrix.")
    return(NULL)
  }
  
  # Using new GSVA 2.0 API (ssgseaParam)
  message("    Preparing GSVA parameters...")
  params <- ssgseaParam(
    exprData     = expr_mat,
    geneSets     = filtered_sigs,
    minSize      = 2,
    maxSize      = 500,
    normalize    = TRUE
  )
  
  message("    Computing scores (this may take a few minutes)...")
  res <- tryCatch(
    gsva(params, verbose = FALSE),
    error = function(e) { 
      message("    GSVA Error: ", e$message)
      # Fallback to old API if new one fails for some reason
      tryCatch(
        gsva(expr_mat, filtered_sigs, method="ssgsea", kcdf="Gaussian", verbose=FALSE),
        error = function(e2) { message("    Fallback Error: ", e2$message); NULL }
      )
    }
  )
  
  if (!is.null(res)) message("    ssGSEA complete. Result dimensions: ", nrow(res), "x", ncol(res))
  res
}

# ── Analysis: Correlate Immune with Risk
correlate_immune_risk <- function(ssgsea_res, risk_df) {
  if (is.null(ssgsea_res) || is.null(risk_df)) return(NULL)
  
  # Match samples
  common_samples <- intersect(colnames(ssgsea_res), risk_df$sample_id)
  if (length(common_samples) < 10) return(NULL)
  
  ss_sub   <- t(ssgsea_res[, common_samples])
  risk_sub <- risk_df |> filter(sample_id %in% common_samples) |> arrange(match(sample_id, common_samples))
  
  # Correlation table
  cor_list <- lapply(colnames(ss_sub), function(cell) {
    ct <- cor.test(ss_sub[, cell], risk_sub$risk_score, method = "spearman")
    tibble(
      Cell_Type = cell,
      Spearman_r = ct$estimate,
      P_Value = ct$p.value
    )
  }) |> bind_rows() |> arrange(Spearman_r)
  
  list(cor_table = cor_list, ss_data = ss_sub, risk_data = risk_sub)
}

# ── Plot: Correlation Barplot
save_immune_risk_barplot <- function(cor_df, cancer, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  
  p <- ggplot(cor_df, aes(x = reorder(Cell_Type, Spearman_r), y = Spearman_r, fill = Spearman_r)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient2(low = "#2171b5", mid = "white", high = "#de2d26", midpoint = 0) +
    labs(
      title = paste(tools::toTitleCase(cancer), "— Immune vs. Risk Score"),
      subtitle = "Spearman Correlation (ssGSEA vs. LASSO-Cox Risk)",
      x = NULL, y = "Spearman r"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
  ggsave(path, plot = p, width = 7, height = 8, dpi = 600, bg = "white")
}

# ── Plot: Heatmap of Immune Scores
save_immune_heatmap <- function(ss_mat, risk_df, cancer, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  
  # Annotate by Risk Group
  ann_col <- risk_df |> 
    select(sample_id, risk_group) |> 
    column_to_rownames("sample_id")
  
  # Z-score for heatmap
  mat <- t(scale(t(ss_mat)))
  
  pheatmap(
    mat,
    annotation_col = ann_col,
    show_colnames = FALSE,
    scale = "none",
    clustering_method = "ward.D2",
    color = colorRampPalette(c("#2171b5", "white", "#de2d26"))(100),
    main = paste(tools::toTitleCase(cancer), "— Immune Infiltration Landscape"),
    annotation_colors = list(risk_group = c(Low = "#3182bd", High = "#de2d26")),
    filename = path,
    width = 10, height = 8
  )
}
