library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
library(tidyverse)
library(rio)

# AnnotationDbi (loaded via org.Hs.eg.db) masks dplyr::select — pin it back
select <- dplyr::select
filter <- dplyr::filter

CANCERS <- c("colorectal", "esophagus", "kidney", "liver", "pancreatic")

# ── Thresholds
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1      # |pooled_log2FC| threshold for ORA

# ── Load meta-analysis results for one cancer
load_meta <- function(cancer) {
  path <- paste0("outputs/", cancer, "/meta_analysis/", cancer, "_meta_analysis_results.rds")
  readRDS(path)
}


# ── Symbol → Entrez conversion
symbols_to_entrez <- function(symbols) {
  bitr(
    geneID  = symbols,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
}

# ── Build ranked list for GSEA
# Rank metric: sign(LFC) * -log10(rem_pvalue)  — captures direction + significance
build_ranked_list <- function(meta_df) {
  df <- meta_df |>
    filter(!is.na(rem_pvalue), !is.na(pooled_log2FC), rem_pvalue > 0) |>
    mutate(rank_score = sign(pooled_log2FC) * -log10(rem_pvalue))
  entrez_map <- symbols_to_entrez(df$Symbol)
  df <- df |>
    left_join(entrez_map, by = c("Symbol" = "SYMBOL")) |>
    filter(!is.na(ENTREZID), !duplicated(ENTREZID)) |>
    arrange(desc(rank_score))
  setNames(df$rank_score, df$ENTREZID)
}

# ── ORA: GO (BP/MF/CC) + KEGG + Reactome
run_ora <- function(entrez_ids, universe_entrez) {

  if (length(entrez_ids) < 5) {
    message("    Too few genes for ORA (n=", length(entrez_ids), "), skipping.")
    return(NULL)
  }

  results <- list()

  for (ont in c("BP", "MF", "CC")) {
    res <- tryCatch(
      enrichGO(
        gene          = entrez_ids,
        universe      = universe_entrez,
        OrgDb         = org.Hs.eg.db,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      ),
      error = function(e) { message("    GO-", ont, " error: ", e$message); NULL }
    )
    results[[paste0("GO_", ont)]] <- res
  }

  results[["KEGG"]] <- tryCatch(
    enrichKEGG(
      gene          = entrez_ids,
      universe      = universe_entrez,
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    ),
    error = function(e) { message("    KEGG error: ", e$message); NULL }
  )

  results[["Reactome"]] <- tryCatch(
    enrichPathway(
      gene          = entrez_ids,
      universe      = universe_entrez,
      organism      = "human",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      readable      = TRUE
    ),
    error = function(e) { message("    Reactome error: ", e$message); NULL }
  )

  results
}

# ── GSEA: GO-BP + KEGG + Hallmark
run_gsea <- function(ranked_list) {

  if (length(ranked_list) < 100) {
    message("    Too few ranked genes (n=", length(ranked_list), "), skipping GSEA.")
    return(NULL)
  }

  results <- list()

  results[["GSEA_GO_BP"]] <- tryCatch(
    gseGO(
      geneList     = ranked_list,
      OrgDb        = org.Hs.eg.db,
      ont          = "BP",
      minGSSize    = 10,
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      verbose      = FALSE
    ),
    error = function(e) { message("    GSEA GO error: ", e$message); NULL }
  )

  results[["GSEA_KEGG"]] <- tryCatch(
    gseKEGG(
      geneList     = ranked_list,
      organism     = "hsa",
      minGSSize    = 10,
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      verbose      = FALSE
    ),
    error = function(e) { message("    GSEA KEGG error: ", e$message); NULL }
  )

  hallmark <- msigdbr(species = "Homo sapiens", category = "H") |>
    select(gs_name, entrez_gene)

  results[["GSEA_Hallmark"]] <- tryCatch(
    GSEA(
      geneList     = ranked_list,
      TERM2GENE    = hallmark,
      minGSSize    = 10,
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      verbose      = FALSE
    ),
    error = function(e) { message("    GSEA Hallmark error: ", e$message); NULL }
  )

  results
}

# ── Export enrichment results to CSV
export_enrichment <- function(result_list, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (name in names(result_list)) {
    res <- result_list[[name]]
    if (is.null(res)) next
    df <- tryCatch(as.data.frame(res), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    export(df, file.path(out_dir, paste0(prefix, "_", name, ".csv")))
  }
}

# ── Plot: dotplot
save_dotplot <- function(res, title, path, width = 10, height = 8) {
  if (is.null(res)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(res), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  p <- dotplot(res, showCategory = 20) + ggtitle(title)
  ggsave(path, plot = p, width = width, height = height, dpi = 300)
  invisible(p)
}

# ── Plot: GSEA enrichment plot (top 3 pathways)
save_gsea_plot <- function(res, title, path, n = 3) {
  if (is.null(res)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(res), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  ids <- head(df$ID, n)
  p <- gseaplot2(res, geneSetID = ids, title = title)
  ggsave(path, plot = p, width = 12, height = 6 * ceiling(n / 2), dpi = 300)
  invisible(p)
}
