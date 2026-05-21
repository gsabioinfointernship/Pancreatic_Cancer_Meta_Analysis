library(STRINGdb)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidyverse)
library(rio)
library(scales)

# AnnotationDbi masks dplyr::select — pin it back
select <- dplyr::select
filter <- dplyr::filter

CANCERS <- c("colorectal", "esophagus", "kidney", "liver", "pancreatic")

# ── Thresholds
PADJ_CUTOFF      <- 0.05
LFC_CUTOFF       <- 1       # |pooled_log2FC| for ORA / PPI node inclusion
STRING_SCORE     <- 700     # combined score threshold (0-1000)
STRING_VERSION   <- "12.0"
HUB_QUANTILE     <- 0.90    # top 10% degree = hub gene
MIN_MODULE_SIZE  <- 5       # skip tiny modules in subgraph plots
MIN_CANCERS_PC   <- 3       # pan-cancer: gene must be sig in >= N cancers

# STRINGdb object — STRINGdb caches downloaded files to working directory
get_string_db <- function() {
  STRINGdb$new(
    version         = STRING_VERSION,
    species         = 9606,
    score_threshold = STRING_SCORE
  )
}

# ── Load meta-analysis results for one cancer
load_meta <- function(cancer) {
  path <- paste0("outputs/", cancer, "/meta_analysis/", cancer, "_meta_analysis_results.rds")
  readRDS(path)
}

# ── Map gene symbols to STRING IDs
map_to_string <- function(gene_df, string_db) {
  # gene_df must have columns: Symbol, pooled_log2FC
  mapped <- string_db$map(
    as.data.frame(gene_df),
    my_data_frame_id_col_names = "Symbol",
    removeUnmappedRows = TRUE
  )
  mapped
}

# ── Build igraph PPI network from STRING interactions
build_ppi_graph <- function(mapped_df, string_db) {
  if (nrow(mapped_df) < 2) return(NULL)

  interactions <- tryCatch(
    string_db$get_interactions(mapped_df$STRING_id),
    error = function(e) { message("  STRING interactions error: ", e$message); NULL }
  )
  if (is.null(interactions) || nrow(interactions) == 0) return(NULL)

  # Resolve STRING IDs back to gene symbols
  id_to_symbol <- setNames(mapped_df$Symbol, mapped_df$STRING_id)
  id_to_lfc    <- setNames(mapped_df$pooled_log2FC, mapped_df$STRING_id)

  edges <- interactions |>
    mutate(
      gene_from = id_to_symbol[from],
      gene_to   = id_to_symbol[to],
      lfc_from  = id_to_lfc[from],
      lfc_to    = id_to_lfc[to]
    ) |>
    filter(!is.na(gene_from), !is.na(gene_to))

  if (nrow(edges) == 0) return(NULL)

  vertex_df <- mapped_df |>
    select(Symbol, pooled_log2FC) |>
    distinct(Symbol, .keep_all = TRUE)

  g <- graph_from_data_frame(
    d        = edges |> select(gene_from, gene_to, combined_score),
    vertices = vertex_df,
    directed = FALSE
  )
  g <- simplify(g)
  g
}

# ── Compute network centrality metrics
compute_metrics <- function(g) {
  tibble(
    Symbol      = V(g)$name,
    Degree      = degree(g),
    Betweenness = betweenness(g, normalized = TRUE),
    Closeness   = closeness(g, normalized = TRUE),
    Eigenvector = eigen_centrality(g)$vector,
    Log2FC      = V(g)$pooled_log2FC
  ) |>
    arrange(desc(Degree))
}

# ── Identify hub genes (top N% by degree)
get_hub_genes <- function(metrics_df, quantile_cut = HUB_QUANTILE) {
  thresh <- quantile(metrics_df$Degree, quantile_cut)
  metrics_df |> filter(Degree >= thresh)
}

# ── Detect network communities (Louvain)
detect_communities <- function(g) {
  comm <- cluster_louvain(g)
  V(g)$community <- membership(comm)
  list(graph = g, communities = comm)
}

# ── Set visual attributes
set_visual_attrs <- function(g, metrics_df) {
  hub_thresh <- quantile(metrics_df$Degree, HUB_QUANTILE)
  V(g)$size  <- rescale(degree(g), to = c(3, 15))
  V(g)$color <- ifelse(V(g)$pooled_log2FC > 0, "#de2d26", "#2171b5")
  V(g)$label <- ifelse(degree(g) >= hub_thresh, V(g)$name, NA)
  g
}

# ── Save full PPI network plot
save_network_plot <- function(g, title, path, seed = 42) {
  if (is.null(g) || vcount(g) < 2) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  png(path, width = 12, height = 12, units = "in", res = 300)
  set.seed(seed)
  plot(
    g,
    vertex.size        = V(g)$size,
    vertex.color       = V(g)$color,
    vertex.label       = V(g)$label,
    vertex.label.cex   = 0.65,
    vertex.label.color = "black",
    edge.width         = 0.4,
    edge.color         = "gray75",
    layout             = layout_with_fr(g),
    main               = title
  )
  legend("topright",
         legend = c("Upregulated", "Downregulated"),
         col    = c("#de2d26", "#2171b5"),
         pch    = 19, bty = "n", pt.cex = 1.5)
  dev.off()
  invisible(NULL)
}

# ── Save ggraph network plot (publication quality)
save_ggraph_plot <- function(g, metrics_df, title, path, seed = 42) {
  if (is.null(g) || vcount(g) < 2) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  hub_thresh <- quantile(metrics_df$Degree, HUB_QUANTILE)
  node_df <- metrics_df |>
    mutate(
      is_hub     = Degree >= hub_thresh,
      Direction  = ifelse(Log2FC > 0, "Up", "Down"),
      label_name = ifelse(is_hub, Symbol, NA_character_)
    )

  set.seed(seed)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.25, color = "gray60", width = 0.3) +
    geom_node_point(
      aes(size = node_df$Degree, color = node_df$Direction),
      alpha = 0.85
    ) +
    geom_node_text(
      aes(label = node_df$label_name),
      repel = TRUE, size = 2.5, max.overlaps = 20, na.rm = TRUE
    ) +
    scale_color_manual(
      values = c(Up = "#de2d26", Down = "#2171b5"),
      name   = "Direction"
    ) +
    scale_size_continuous(range = c(1.5, 8), name = "Degree") +
    labs(title = title) +
    theme_graph(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = "right"
    )

  ggsave(path, plot = p, width = 14, height = 12, dpi = 300)
  invisible(p)
}

# ── Hub gene barplot
save_hub_barplot <- function(hub_df, title, path, n = 25) {
  if (is.null(hub_df) || nrow(hub_df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- head(hub_df, n)
  p <- ggplot(df, aes(x = reorder(Symbol, Degree), y = Degree, fill = Log2FC)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(
      low      = "#2171b5",
      mid      = "white",
      high     = "#de2d26",
      midpoint = 0,
      name     = "Log2FC"
    ) +
    coord_flip() +
    labs(title = title, x = NULL, y = "Degree Centrality") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(path, plot = p, width = 8, height = 7, dpi = 300)
  invisible(p)
}

# ── Module subgraph plots (top N modules)
save_module_plots <- function(g, fig_dir, prefix, n_modules = 5, seed = 42) {
  if (is.null(g) || is.null(V(g)$community)) return(invisible(NULL))
  module_sizes <- sort(table(V(g)$community), decreasing = TRUE)
  top_ids      <- names(head(module_sizes, n_modules))

  for (mid in top_ids) {
    vids <- which(V(g)$community == mid)
    subg <- induced_subgraph(g, vids = vids)
    if (vcount(subg) < MIN_MODULE_SIZE) next

    path <- file.path(fig_dir, paste0(prefix, "_module_", mid, "_network.png"))
    png(path, width = 8, height = 8, units = "in", res = 300)
    set.seed(seed)
    sub_deg <- degree(subg)
    plot(
      subg,
      vertex.size        = rescale(sub_deg, to = c(5, 18)),
      vertex.color       = V(subg)$color,
      vertex.label       = V(subg)$name,
      vertex.label.cex   = 0.6,
      vertex.label.color = "black",
      edge.width         = 1,
      edge.color         = "gray70",
      layout             = layout_with_fr(subg),
      main               = paste0(prefix, " — Module ", mid,
                                  " (", vcount(subg), " genes)")
    )
    dev.off()
  }
}

# ── Save metrics + hub genes to CSV
export_network <- function(metrics_df, hub_df, modules_df, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  export(metrics_df, file.path(out_dir, paste0(prefix, "_network_metrics.csv")))
  export(hub_df,     file.path(out_dir, paste0(prefix, "_hub_genes.csv")))
  export(modules_df, file.path(out_dir, paste0(prefix, "_network_modules.csv")))
}
