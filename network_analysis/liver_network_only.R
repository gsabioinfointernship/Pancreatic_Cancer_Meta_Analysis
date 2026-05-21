# ============================================================================
# Liver-only network analysis (30 DEGs — targeted run)
# Lowered STRING_SCORE to 400 for sparse gene set to capture more interactions
# ============================================================================

source("network_analysis/00_network_functions.R")

cancer <- "liver"

out_csv <- paste0("outputs/", cancer, "/network/")
out_fig <- paste0("outputs/", cancer, "/figures/network/")
dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

message("=== LIVER NETWORK ANALYSIS (targeted) ===")

# Load meta results
meta <- load_meta(cancer)
sig <- meta |>
  filter(rem_padj < PADJ_CUTOFF, abs(pooled_log2FC) >= LFC_CUTOFF) |>
  select(Symbol, pooled_log2FC) |>
  filter(!is.na(Symbol), Symbol != "")

message("Significant DEGs: ", nrow(sig))

# Try STRING at primary threshold (700) first, fall back to 400 if needed
for (score_thresh in c(700, 400)) {
  message("\nTrying STRING score >= ", score_thresh, " ...")

  string_db_local <- STRINGdb$new(
    version         = STRING_VERSION,
    species         = 9606,
    score_threshold = score_thresh
  )

  mapped <- tryCatch(
    map_to_string(sig, string_db_local),
    error = function(e) { message("  Mapping error: ", e$message); NULL }
  )

  if (is.null(mapped) || nrow(mapped) < 2) {
    message("  Insufficient mappings at ", score_thresh)
    next
  }
  message("  Mapped: ", nrow(mapped), " genes")

  g <- build_ppi_graph(mapped, string_db_local)
  if (!is.null(g) && vcount(g) >= 3) {
    message("  Network: ", vcount(g), " nodes, ", ecount(g), " edges")
    break
  } else {
    message("  No usable graph at ", score_thresh)
    g <- NULL
  }
}

if (is.null(g)) {
  message("\nNo STRING interactions found for liver DEGs.")
  message("Writing empty hub_genes.csv as placeholder.")

  empty_hub <- data.frame(
    Symbol = character(0), Degree = integer(0),
    Betweenness = numeric(0), Closeness = numeric(0),
    Eigenvector = numeric(0), Log2FC = numeric(0)
  )
  write.csv(empty_hub, file.path(out_csv, paste0(cancer, "_hub_genes.csv")),
            row.names = FALSE)

  empty_metrics <- data.frame(
    Symbol = sig$Symbol,
    Degree = 0L, Betweenness = 0, Closeness = 0,
    Eigenvector = 0, Log2FC = sig$pooled_log2FC
  )
  write.csv(empty_metrics,
            file.path(out_csv, paste0(cancer, "_network_metrics.csv")),
            row.names = FALSE)

  message("Done — liver network is genuinely sparse (30 DEGs × STRING >= 400).")
  message("Recommend: acknowledge as study limitation in manuscript.")
  quit(save = "no", status = 0)
}

# Centrality + hub genes
metrics   <- compute_metrics(g)
hub_genes <- get_hub_genes(metrics)
message("Hub genes (top 10%): ", nrow(hub_genes))

# Community detection
comm_result <- detect_communities(g)
g <- comm_result$graph

modules_df <- tibble(
  Symbol  = V(g)$name,
  Module  = V(g)$community,
  Degree  = degree(g),
  Log2FC  = V(g)$pooled_log2FC
) |> arrange(Module, desc(Degree))

# Export
export_network(metrics, hub_genes, modules_df, out_csv, cancer)
message("CSVs saved to: ", out_csv)

# Visual attributes
g <- set_visual_attrs(g, metrics)

# Figures
save_network_plot(g,
  title = "Liver PPI Network",
  path  = file.path(out_fig, paste0(cancer, "_PPI_network.png"))
)

save_ggraph_plot(g, metrics,
  title = "Liver PPI Network",
  path  = file.path(out_fig, paste0(cancer, "_PPI_ggraph.png"))
)

save_hub_barplot(hub_genes,
  title = "Liver — Top Hub Genes",
  path  = file.path(out_fig, paste0(cancer, "_hub_genes_barplot.png"))
)

save_module_plots(g, out_fig, cancer)

message("\n=== LIVER NETWORK COMPLETE ===")
message("Nodes: ", vcount(g), "  Edges: ", ecount(g),
        "  Hub genes: ", nrow(hub_genes))
