source("network_analysis/00_network_functions.R")

# Initialise STRINGdb once (shared cache across cancers)
string_db <- get_string_db()

for (cancer in CANCERS) {

  message("\n=== ", toupper(cancer), " ===")

  # ── Directories
  out_csv <- paste0("outputs/", cancer, "/network/")
  out_fig <- paste0("outputs/", cancer, "/figures/network/")
  dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

  # ── Load meta results
  meta <- tryCatch(load_meta(cancer), error = function(e) {
    message("  SKIPPING — could not load meta: ", e$message); NULL
  })
  if (is.null(meta)) next

  # ── Significant DEGs
  sig <- meta |>
    filter(rem_padj < PADJ_CUTOFF, abs(pooled_log2FC) >= LFC_CUTOFF) |>
    select(Symbol, pooled_log2FC) |>
    filter(!is.na(Symbol), Symbol != "")

  message("  Significant DEGs: ", nrow(sig))

  if (nrow(sig) < 5) {
    message("  SKIPPING — too few significant DEGs")
    next
  }

  # ── Map to STRING
  message("  Mapping to STRING...")
  mapped <- tryCatch(
    map_to_string(sig, string_db),
    error = function(e) { message("  STRING mapping error: ", e$message); NULL }
  )
  if (is.null(mapped) || nrow(mapped) < 2) {
    message("  SKIPPING — insufficient STRING mappings")
    next
  }
  message("  Mapped genes: ", nrow(mapped))

  # ── Build PPI graph
  message("  Building PPI network...")
  g <- build_ppi_graph(mapped, string_db)
  if (is.null(g)) {
    message("  SKIPPING — no interactions found")
    next
  }
  message("  Nodes: ", vcount(g), "  Edges: ", ecount(g))

  # ── Centrality metrics
  metrics <- compute_metrics(g)
  hub_genes <- get_hub_genes(metrics)
  message("  Hub genes (top ", round((1 - HUB_QUANTILE) * 100), "%): ", nrow(hub_genes))

  # ── Community detection
  comm_result <- detect_communities(g)
  g <- comm_result$graph

  modules_df <- tibble(
    Symbol    = V(g)$name,
    Module    = V(g)$community,
    Degree    = degree(g),
    Log2FC    = V(g)$pooled_log2FC
  ) |> arrange(Module, desc(Degree))

  # ── Export CSVs
  export_network(metrics, hub_genes, modules_df, out_csv, cancer)

  # ── Visual attributes
  g <- set_visual_attrs(g, metrics)

  # ── Figures
  message("  Saving figures...")

  # Base PPI network
  save_network_plot(
    g     = g,
    title = paste(tools::toTitleCase(cancer), "PPI Network"),
    path  = file.path(out_fig, paste0(cancer, "_PPI_network.png"))
  )

  # ggraph (publication quality)
  save_ggraph_plot(
    g          = g,
    metrics_df = metrics,
    title      = paste(tools::toTitleCase(cancer), "PPI Network"),
    path       = file.path(out_fig, paste0(cancer, "_PPI_ggraph.png"))
  )

  # Hub gene barplot
  save_hub_barplot(
    hub_df = hub_genes,
    title  = paste(tools::toTitleCase(cancer), "— Top Hub Genes"),
    path   = file.path(out_fig, paste0(cancer, "_hub_genes_barplot.png"))
  )

  # Module subgraphs
  save_module_plots(
    g       = g,
    fig_dir = out_fig,
    prefix  = cancer
  )

  message("  Done: ", cancer)
}

message("\nPer-cancer network analysis complete.")
