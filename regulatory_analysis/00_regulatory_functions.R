library(tidyverse)
library(clusterProfiler)
library(rio)
library(ggplot2)
library(ggrepel)

# ── Shared Configuration
CANCERS         <- c("colorectal", "esophagus", "kidney", "liver", "pancreatic")
PADJ_CUTOFF     <- 0.05
LFC_CUTOFF      <- 1

# ── Helper: Load Meta-Analysis Results
load_meta <- function(cancer) {
  path <- paste0("outputs/", cancer, "/meta_analysis/", cancer, "_meta_analysis_results.rds")
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

# ── Core: Master Regulator Analysis (ORA for TFs)
# Using TRRUST v2 logic: Enrichment of sig DEGs in known TF-target sets
run_tf_ora <- function(meta_df, tf_db, direction = "up") {
  
  if (direction == "up") {
    sig_genes <- meta_df |> filter(rem_padj < PADJ_CUTOFF, pooled_log2FC >= LFC_CUTOFF) |> pull(Symbol)
  } else {
    sig_genes <- meta_df |> filter(rem_padj < PADJ_CUTOFF, pooled_log2FC <= -LFC_CUTOFF) |> pull(Symbol)
  }
  
  if (length(sig_genes) < 5) return(NULL)
  
  # ORA using clusterProfiler::enricher
  # tf_db should be a long format df: TF, Target
  res <- enricher(
    gene         = sig_genes,
    TERM2GENE    = tf_db,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize    = 3,
    maxGSSize    = 1000
  )
  
  if (is.null(res)) return(NULL)
  as.data.frame(res)
}

# ── Core: Master Regulator Analysis (GSEA for TFs)
run_tf_gsea <- function(meta_df, tf_db) {
  # Rank by sign(LFC) * -log10(p)
  ranked <- meta_df |> 
    mutate(score = sign(pooled_log2FC) * -log10(rem_pvalue)) |>
    filter(!is.na(score), !is.infinite(score)) |>
    arrange(desc(score))
    
  gene_list <- ranked$score
  names(gene_list) <- ranked$Symbol
  
  res <- GSEA(
    geneList     = gene_list,
    TERM2GENE    = tf_db,
    pvalueCutoff = 0.05,
    minGSSize    = 3,
    maxGSSize    = 500,
    verbose      = FALSE
  )
  
  if (is.null(res)) return(NULL)
  as.data.frame(res)
}

# ── Plot: TF Enrichment Dotplot
save_tf_dotplot <- function(tf_res, cancer, direction, path) {
  if (is.null(tf_res) || nrow(tf_res) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  
  plot_df <- tf_res |> 
    head(15) |>
    mutate(ID = factor(ID, levels = rev(ID)))
    
  p <- ggplot(plot_df, aes(x = Count, y = ID, size = -log10(p.adjust), color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#de2d26", high = "#fee0d2") +
    labs(
      title = paste(tools::toTitleCase(cancer), "— Master Regulators"),
      subtitle = paste("TF Enrichment among", direction, "regulated DEGs"),
      x = "Target Gene Count", y = "Transcription Factor",
      size = "-log10(FDR)", color = "Adj. P-value"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
  ggsave(path, plot = p, width = 7, height = 6, dpi = 600, bg = "white")
}
