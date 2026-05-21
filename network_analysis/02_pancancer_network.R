source("network_analysis/00_network_functions.R")

message("=== PAN-CANCER NETWORK ANALYSIS ===")

# ── Directories
out_csv <- "outputs/pancancer/network/"
out_fig <- "outputs/pancancer/figures/network/"
dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

# Initialise STRINGdb (shared cache)
string_db <- get_string_db()

# ============================================================================
# 1. Load all per-cancer meta results
# ============================================================================

all_meta <- lapply(setNames(CANCERS, CANCERS), function(cancer) {
  tryCatch(load_meta(cancer), error = function(e) {
    message("  Could not load ", cancer, ": ", e$message); NULL
  })
})
all_meta <- Filter(Negate(is.null), all_meta)
message("Loaded meta for: ", paste(names(all_meta), collapse = ", "))

# ============================================================================
# 2. Pan-cancer consistent DEGs
# ============================================================================

sig_per_cancer <- lapply(all_meta, function(meta) {
  meta |>
    filter(rem_padj < PADJ_CUTOFF, abs(pooled_log2FC) >= LFC_CUTOFF) |>
    select(Symbol, pooled_log2FC)
})

gene_counts <- bind_rows(sig_per_cancer, .id = "cancer") |>
  group_by(Symbol) |>
  summarise(
    n_cancers  = n(),
    n_up       = sum(pooled_log2FC > 0),
    n_down     = sum(pooled_log2FC < 0),
    mean_lfc   = mean(pooled_log2FC),
    sd_lfc     = sd(pooled_log2FC),
    cancers    = paste(cancer, collapse = "|"),
    .groups    = "drop"
  ) |>
  filter(n_cancers >= MIN_CANCERS_PC) |>
  arrange(desc(n_cancers), desc(abs(mean_lfc)))

message("Pan-cancer consistent DEGs (>=", MIN_CANCERS_PC, " cancers): ", nrow(gene_counts))

pancancer_up   <- gene_counts |> filter(n_down == 0) |> pull(Symbol)
pancancer_down <- gene_counts |> filter(n_up   == 0) |> pull(Symbol)
pancancer_all  <- gene_counts$Symbol

message("  Consistent up: ",   length(pancancer_up),
        "  Consistent down: ", length(pancancer_down))

export(gene_counts, file.path(out_csv, "pancancer_consistent_DEGs.csv"))

# ============================================================================
# 3. Pan-cancer PPI network (consistent DEGs)
# ============================================================================

if (length(pancancer_all) >= 5) {
  message("Building pan-cancer PPI network...")

  pc_gene_df <- gene_counts |>
    select(Symbol, mean_lfc) |>
    rename(pooled_log2FC = mean_lfc)

  pc_mapped <- tryCatch(
    map_to_string(pc_gene_df, string_db),
    error = function(e) { message("  STRING mapping error: ", e$message); NULL }
  )

  if (!is.null(pc_mapped) && nrow(pc_mapped) >= 2) {
    message("  Mapped: ", nrow(pc_mapped), " genes")

    g_pc <- build_ppi_graph(pc_mapped, string_db)

    if (!is.null(g_pc)) {
      message("  Nodes: ", vcount(g_pc), "  Edges: ", ecount(g_pc))

      metrics_pc   <- compute_metrics(g_pc)
      hub_genes_pc <- get_hub_genes(metrics_pc)

      comm_pc <- detect_communities(g_pc)
      g_pc    <- comm_pc$graph

      # Add cancer membership info to vertex attributes
      cancer_membership <- gene_counts |>
        select(Symbol, n_cancers, cancers)
      V(g_pc)$n_cancers <- cancer_membership$n_cancers[
        match(V(g_pc)$name, cancer_membership$Symbol)
      ]
      V(g_pc)$cancers_list <- cancer_membership$cancers[
        match(V(g_pc)$name, cancer_membership$Symbol)
      ]

      modules_pc <- tibble(
        Symbol    = V(g_pc)$name,
        Module    = V(g_pc)$community,
        N_Cancers = V(g_pc)$n_cancers,
        Degree    = degree(g_pc),
        Log2FC    = V(g_pc)$pooled_log2FC
      ) |> arrange(Module, desc(Degree))

      # Add n_cancers to metrics
      metrics_pc <- metrics_pc |>
        left_join(gene_counts |> select(Symbol, n_cancers, cancers), by = "Symbol")

      export_network(metrics_pc, hub_genes_pc, modules_pc, out_csv, "pancancer")

      message("  Hub genes (top ",
              round((1 - HUB_QUANTILE) * 100), "%): ", nrow(hub_genes_pc))

      # Visual attrs
      g_pc <- set_visual_attrs(g_pc, metrics_pc)

      # Full PPI network
      save_network_plot(
        g     = g_pc,
        title = "Pan-Cancer PPI Network (Consistent DEGs)",
        path  = file.path(out_fig, "pancancer_PPI_network.png")
      )

      # ggraph — nodes sized by n_cancers
      hub_thresh <- quantile(metrics_pc$Degree, HUB_QUANTILE)
      node_df <- metrics_pc |>
        mutate(
          Direction  = ifelse(Log2FC > 0, "Up", "Down"),
          label_name = ifelse(Degree >= hub_thresh, Symbol, NA_character_)
        )

      set.seed(42)
      p_pc <- ggraph(g_pc, layout = "fr") +
        geom_edge_link(alpha = 0.2, color = "gray60", width = 0.3) +
        geom_node_point(
          aes(size = node_df$n_cancers, color = node_df$Direction),
          alpha = 0.85
        ) +
        geom_node_text(
          aes(label = node_df$label_name),
          repel = TRUE, size = 2.5, max.overlaps = 25, na.rm = TRUE
        ) +
        scale_color_manual(
          values = c(Up = "#de2d26", Down = "#2171b5"),
          name   = "Direction"
        ) +
        scale_size_continuous(
          range  = c(2, 10),
          name   = "# Cancers",
          breaks = seq(MIN_CANCERS_PC, length(all_meta))
        ) +
        labs(title = "Pan-Cancer PPI Network — Consistent DEGs") +
        theme_graph(base_family = "sans") +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
          legend.position = "right"
        )
      ggsave(file.path(out_fig, "pancancer_PPI_ggraph.png"),
             plot = p_pc, width = 14, height = 12, dpi = 300)

      # Hub gene barplot
      save_hub_barplot(
        hub_df = hub_genes_pc,
        title  = "Pan-Cancer Hub Genes",
        path   = file.path(out_fig, "pancancer_hub_genes_barplot.png")
      )

      # Module subgraphs
      save_module_plots(g_pc, out_fig, "pancancer")
    }
  }
} else {
  message("Too few pan-cancer consistent DEGs for network (n=", length(pancancer_all), ")")
}

# ============================================================================
# 4. Cross-cancer hub gene comparison
# ============================================================================

message("Building cross-cancer hub gene comparison...")

# Load per-cancer hub gene files (generated by 01_per_cancer_network.R)
hub_per_cancer <- lapply(setNames(CANCERS, CANCERS), function(cancer) {
  path <- paste0("outputs/", cancer, "/network/", cancer, "_hub_genes.csv")
  if (!file.exists(path)) return(NULL)
  import(path) |> mutate(cancer = cancer)
})
hub_per_cancer <- Filter(Negate(is.null), hub_per_cancer)

if (length(hub_per_cancer) >= 2) {
  all_hubs <- bind_rows(hub_per_cancer)

  # Genes that are hubs in >= 2 cancers
  shared_hubs <- all_hubs |>
    count(Symbol, name = "n_cancers_hub") |>
    filter(n_cancers_hub >= 2) |>
    left_join(
      all_hubs |>
        group_by(Symbol) |>
        summarise(
          mean_degree = mean(Degree),
          mean_lfc    = mean(Log2FC, na.rm = TRUE),
          cancer_list = paste(cancer, collapse = "|"),
          .groups = "drop"
        ),
      by = "Symbol"
    ) |>
    arrange(desc(n_cancers_hub), desc(mean_degree))

  export(shared_hubs, file.path(out_csv, "pancancer_shared_hub_genes.csv"))
  message("  Shared hub genes (>=2 cancers): ", nrow(shared_hubs))

  # Heatmap-style: hub degree per cancer
  hub_matrix <- all_hubs |>
    select(Symbol, cancer, Degree) |>
    pivot_wider(names_from = cancer, values_from = Degree, values_fill = 0)

  export(hub_matrix, file.path(out_csv, "pancancer_hub_degree_matrix.csv"))

  # ── Plot: shared hub genes bubble chart
  if (nrow(shared_hubs) > 0) {
    plot_df <- all_hubs |>
      filter(Symbol %in% shared_hubs$Symbol) |>
      mutate(Direction = ifelse(Log2FC > 0, "Up", "Down"))

    p_shared <- ggplot(plot_df, aes(x = cancer, y = reorder(Symbol, Degree),
                                    size = Degree, color = Direction)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c(Up = "#de2d26", Down = "#2171b5")) +
      scale_size_continuous(range = c(2, 10), name = "Degree") +
      labs(
        title = "Shared Hub Genes Across Cancers",
        x     = "Cancer Type",
        y     = "Gene Symbol"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x  = element_text(angle = 35, hjust = 1),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_line(color = "gray90")
      )
    ggsave(file.path(out_fig, "pancancer_shared_hubs_bubble.png"),
           plot = p_shared, width = 10,
           height = min(48, max(6, nrow(shared_hubs) * 0.3 + 3)),
           dpi = 300)
  }

  # ── Plot: hub gene count per cancer
  p_count <- all_hubs |>
    count(cancer, name = "n_hubs") |>
    ggplot(aes(x = reorder(cancer, n_hubs), y = n_hubs, fill = cancer)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_text(aes(label = n_hubs), hjust = -0.2, size = 3.5) +
    coord_flip() +
    labs(title = "Hub Gene Count per Cancer", x = NULL, y = "# Hub Genes") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(out_fig, "pancancer_hub_count_per_cancer.png"),
         plot = p_count, width = 7, height = 5, dpi = 300)
}

# ============================================================================
# 5. Cross-cancer network overlap: UpSet-style summary
# ============================================================================

message("Building cross-cancer network overlap table...")

# Which cancers each hub gene appears in
hub_presence <- lapply(CANCERS, function(cancer) {
  path <- paste0("outputs/", cancer, "/network/", cancer, "_hub_genes.csv")
  if (!file.exists(path)) return(character(0))
  import(path)$Symbol
})
names(hub_presence) <- CANCERS

all_hub_symbols <- unique(unlist(hub_presence))

if (length(all_hub_symbols) > 0) {
  presence_df <- tibble(Symbol = all_hub_symbols)
  for (cancer in CANCERS) {
    presence_df[[cancer]] <- as.integer(all_hub_symbols %in% hub_presence[[cancer]])
  }
  presence_df <- presence_df |>
    mutate(total = rowSums(across(all_of(CANCERS)))) |>
    arrange(desc(total))

  export(presence_df, file.path(out_csv, "pancancer_hub_presence_matrix.csv"))
  message("  Hub presence matrix saved: ", nrow(presence_df), " unique hub genes")
}

message("\nPan-cancer network analysis complete.")
