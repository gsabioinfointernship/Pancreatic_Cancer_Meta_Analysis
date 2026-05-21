source("functional_enrichment/00_enrichment_functions.R")

for (cancer in CANCERS) {

  message("\n=== ", toupper(cancer), " ===")

  # ── Directories
  out_csv  <- paste0("outputs/", cancer, "/enrichment/")
  out_fig  <- paste0("figures/", cancer, "/enrichment/")
  dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

  # ── Load meta results
  meta <- load_meta(cancer)
  message("  Genes in meta: ", nrow(meta))

  # ── Universe: all tested genes with Entrez mapping
  universe_map <- symbols_to_entrez(meta$Symbol)
  universe_entrez <- universe_map$ENTREZID

  # ── Significant DEGs for ORA
  sig <- meta |>
    filter(rem_padj < PADJ_CUTOFF, abs(pooled_log2FC) >= LFC_CUTOFF)

  message("  Significant DEGs (padj<", PADJ_CUTOFF, ", |LFC|>=", LFC_CUTOFF, "): ", nrow(sig))

  sig_up   <- filter(sig, pooled_log2FC >  0)
  sig_down <- filter(sig, pooled_log2FC <  0)

  # Convert to Entrez
  entrez_up   <- symbols_to_entrez(sig_up$Symbol)$ENTREZID
  entrez_down <- symbols_to_entrez(sig_down$Symbol)$ENTREZID
  entrez_all  <- symbols_to_entrez(sig$Symbol)$ENTREZID

  message("  Up: ", length(entrez_up), "  Down: ", length(entrez_down))

  # ── ORA
  message("  Running ORA...")

  ora_up   <- run_ora(entrez_up,   universe_entrez)
  ora_down <- run_ora(entrez_down, universe_entrez)
  ora_all  <- run_ora(entrez_all,  universe_entrez)

  export_enrichment(ora_up,   out_csv, paste0(cancer, "_ORA_up"))
  export_enrichment(ora_down, out_csv, paste0(cancer, "_ORA_down"))
  export_enrichment(ora_all,  out_csv, paste0(cancer, "_ORA_all"))

  # Dotplots for all-DEG ORA
  for (name in names(ora_all)) {
    save_dotplot(
      res   = ora_all[[name]],
      title = paste(cancer, name),
      path  = file.path(out_fig, paste0(cancer, "_ORA_all_", name, "_dotplot.png"))
    )
  }

  # ── GSEA
  message("  Running GSEA...")

  ranked <- build_ranked_list(meta)
  message("  Ranked genes for GSEA: ", length(ranked))

  gsea_res <- run_gsea(ranked)

  export_enrichment(gsea_res, out_csv, paste0(cancer, "_GSEA"))

  # GSEA plots
  for (name in names(gsea_res)) {
    save_dotplot(
      res   = gsea_res[[name]],
      title = paste(cancer, name),
      path  = file.path(out_fig, paste0(cancer, "_", name, "_dotplot.png"))
    )
    save_gsea_plot(
      res   = gsea_res[[name]],
      title = paste(cancer, name, "- Top Pathways"),
      path  = file.path(out_fig, paste0(cancer, "_", name, "_enrichplot.png"))
    )
  }

  message("  Done: ", cancer)
}
