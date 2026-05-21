source("functional_enrichment/00_enrichment_functions.R")

message("=== PAN-CANCER ENRICHMENT ===")

out_csv <- "outputs/pancancer/enrichment/"
out_fig <- "figures/pancancer/enrichment/"
dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

# ── Load all meta results
all_meta <- lapply(setNames(CANCERS, CANCERS), load_meta)

# ── Pan-cancer consistent DEGs
# Genes significant in >= 3 cancers, same direction in all where significant
MIN_CANCERS <- 3

sig_per_cancer <- lapply(all_meta, function(meta) {
  meta |>
    filter(rem_padj < PADJ_CUTOFF, abs(pooled_log2FC) >= LFC_CUTOFF) |>
    select(Symbol, pooled_log2FC)
})

# Count occurrences and direction per gene
gene_counts <- bind_rows(sig_per_cancer, .id = "cancer") |>
  group_by(Symbol) |>
  summarise(
    n_cancers  = n(),
    n_up       = sum(pooled_log2FC > 0),
    n_down     = sum(pooled_log2FC < 0),
    mean_lfc   = mean(pooled_log2FC),
    .groups    = "drop"
  ) |>
  filter(n_cancers >= MIN_CANCERS)

message("  Pan-cancer consistent DEGs (>=", MIN_CANCERS, " cancers): ", nrow(gene_counts))

# Direction-consistent: up in all or down in all where significant
pancancer_up   <- gene_counts |> filter(n_down == 0) |> pull(Symbol)
pancancer_down <- gene_counts |> filter(n_up   == 0) |> pull(Symbol)
pancancer_all  <- gene_counts$Symbol

message("  Consistent up: ", length(pancancer_up),
        "  Consistent down: ", length(pancancer_down),
        "  Mixed: ", length(pancancer_all) - length(pancancer_up) - length(pancancer_down))

export(gene_counts, file.path(out_csv, "pancancer_consistent_DEGs.csv"))

# ── Universe: union of all tested genes
all_symbols  <- unique(unlist(lapply(all_meta, function(x) x$Symbol)))
universe_map <- symbols_to_entrez(all_symbols)
universe_entrez <- universe_map$ENTREZID

# ── ORA on pan-cancer gene sets
message("  Running ORA...")

entrez_pc_up   <- symbols_to_entrez(pancancer_up)$ENTREZID
entrez_pc_down <- symbols_to_entrez(pancancer_down)$ENTREZID
entrez_pc_all  <- symbols_to_entrez(pancancer_all)$ENTREZID

ora_up   <- run_ora(entrez_pc_up,   universe_entrez)
ora_down <- run_ora(entrez_pc_down, universe_entrez)
ora_all  <- run_ora(entrez_pc_all,  universe_entrez)

export_enrichment(ora_up,   out_csv, "pancancer_ORA_up")
export_enrichment(ora_down, out_csv, "pancancer_ORA_down")
export_enrichment(ora_all,  out_csv, "pancancer_ORA_all")

for (name in names(ora_all)) {
  save_dotplot(
    res   = ora_all[[name]],
    title = paste("Pan-Cancer", name),
    path  = file.path(out_fig, paste0("pancancer_ORA_all_", name, "_dotplot.png"))
  )
}

# ── Pan-cancer GSEA: meta-rank across cancers
# Aggregate rank score: mean(sign(LFC) * -log10(rem_pvalue)) across cancers
message("  Building pan-cancer ranked list...")

ranked_per_cancer <- lapply(all_meta, function(meta) {
  meta |>
    filter(!is.na(rem_pvalue), !is.na(pooled_log2FC), rem_pvalue > 0) |>
    mutate(rank_score = sign(pooled_log2FC) * -log10(rem_pvalue)) |>
    select(Symbol, rank_score)
})

pancancer_ranks <- bind_rows(ranked_per_cancer) |>
  group_by(Symbol) |>
  summarise(agg_score = mean(rank_score), .groups = "drop")
entrez_map <- symbols_to_entrez(pancancer_ranks$Symbol)
pancancer_ranks <- pancancer_ranks |>
  left_join(entrez_map, by = c("Symbol" = "SYMBOL")) |>
  filter(!is.na(ENTREZID), !duplicated(ENTREZID)) |>
  arrange(desc(agg_score))

ranked_list <- setNames(pancancer_ranks$agg_score, pancancer_ranks$ENTREZID)
message("  Pan-cancer ranked genes: ", length(ranked_list))

gsea_res <- run_gsea(ranked_list)

export_enrichment(gsea_res, out_csv, "pancancer_GSEA")

for (name in names(gsea_res)) {
  save_dotplot(
    res   = gsea_res[[name]],
    title = paste("Pan-Cancer", name),
    path  = file.path(out_fig, paste0("pancancer_", name, "_dotplot.png"))
  )
  save_gsea_plot(
    res   = gsea_res[[name]],
    title = paste("Pan-Cancer", name, "- Top Pathways"),
    path  = file.path(out_fig, paste0("pancancer_", name, "_enrichplot.png"))
  )
}

# ── compareCluster: side-by-side ORA across all 5 cancers
message("  Running compareCluster across cancers...")

gene_list_by_cancer <- lapply(setNames(CANCERS, CANCERS), function(cancer) {
  sig <- all_meta[[cancer]] |>
    filter(rem_padj < PADJ_CUTOFF, abs(pooled_log2FC) >= LFC_CUTOFF) |>
    pull(Symbol)
  entrez <- tryCatch(symbols_to_entrez(sig)$ENTREZID, error = function(e) character(0))
  entrez
})

# Remove cancers with too few genes
gene_list_by_cancer <- Filter(function(x) length(x) >= 5, gene_list_by_cancer)

if (length(gene_list_by_cancer) >= 2) {
  cc_go <- tryCatch(
    compareCluster(
      geneClusters  = gene_list_by_cancer,
      fun           = "enrichGO",
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    ),
    error = function(e) { message("  compareCluster GO error: ", e$message); NULL }
  )

  if (!is.null(cc_go)) {
    export(as.data.frame(cc_go), file.path(out_csv, "compareCluster_GO_BP.csv"))
    p_cc <- dotplot(cc_go, showCategory = 10) +
      ggtitle("GO BP — Cross-Cancer Comparison") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out_fig, "compareCluster_GO_BP_dotplot.png"),
           plot = p_cc, width = 14, height = 10, dpi = 300)
  }

  cc_kegg <- tryCatch(
    compareCluster(
      geneClusters  = gene_list_by_cancer,
      fun           = "enrichKEGG",
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    ),
    error = function(e) { message("  compareCluster KEGG error: ", e$message); NULL }
  )

  if (!is.null(cc_kegg)) {
    export(as.data.frame(cc_kegg), file.path(out_csv, "compareCluster_KEGG.csv"))
    p_cc_kegg <- dotplot(cc_kegg, showCategory = 10) +
      ggtitle("KEGG — Cross-Cancer Comparison") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out_fig, "compareCluster_KEGG_dotplot.png"),
           plot = p_cc_kegg, width = 14, height = 10, dpi = 300)
  }
}
