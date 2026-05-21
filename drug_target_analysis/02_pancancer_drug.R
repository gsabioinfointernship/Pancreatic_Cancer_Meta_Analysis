source("drug_target_analysis/00_drug_functions.R")

message("=== PAN-CANCER DRUG TARGET ANALYSIS ===")

# ── Directories
out_csv <- "outputs/pancancer/drug_targets/"
out_fig <- "outputs/pancancer/figures/drug_targets/"
dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. Load pan-cancer consistent DEGs + shared hub genes
# ============================================================================

pc_degs_path <- "outputs/pancancer/network/pancancer_consistent_DEGs.csv"
shared_hubs_path <- "outputs/pancancer/network/pancancer_shared_hub_genes.csv"

if (!file.exists(pc_degs_path)) {
  stop("Pan-cancer DEGs not found. Run network analysis first: ",
       "network_analysis/02_pancancer_network.R")
}

pc_degs <- import(pc_degs_path)
message("Pan-cancer consistent DEGs: ", nrow(pc_degs))

shared_hubs <- if (file.exists(shared_hubs_path)) {
  import(shared_hubs_path)
} else {
  message("Shared hub file not found — hub flagging disabled")
  NULL
}

# Build a minimal meta-like data frame for pan-cancer DEGs
pc_meta <- pc_degs |>
  rename(pooled_log2FC = mean_lfc) |>
  mutate(rem_pvalue = NA_real_, rem_padj = 0)  # padj=0 → all pass filter

# ============================================================================
# 2. DGIdb query on pan-cancer consistent DEGs
# ============================================================================

message("Querying DGIdb for pan-cancer consistent DEGs...")
iact_pc <- query_dgidb(pc_degs$Symbol)

if (!is.null(iact_pc) && nrow(iact_pc) > 0) {
  message("Drug interactions found: ", nrow(iact_pc))
  message("Druggable pan-cancer genes: ", n_distinct(iact_pc$Symbol))

  hub_stub <- if (!is.null(shared_hubs)) {
    shared_hubs |> select(Symbol) |> rename(Symbol = Symbol)
  } else NULL

  priority_pc <- build_priority_table(iact_pc, pc_meta, hub_stub)

  # Add pan-cancer specific columns
  priority_pc <- priority_pc |>
    left_join(pc_degs |> select(Symbol, n_cancers, cancers), by = "Symbol")

  export_drug(iact_pc, priority_pc, out_csv, "pancancer")

  high_priority_pc <- filter(priority_pc, Priority == "High")
  message("High-priority pan-cancer targets: ", nrow(high_priority_pc))
  if (nrow(high_priority_pc) > 0)
    safe_export(high_priority_pc,
                file.path(out_csv, "pancancer_high_priority_targets.csv"))

  # ── Figures
  save_drug_lollipop(
    priority_df = priority_pc,
    title       = "Pan-Cancer — Top Druggable Targets (Consistent DEGs)",
    path        = file.path(out_fig, "pancancer_drug_lollipop.png")
  )

  save_family_barplot(
    priority_df = priority_pc,
    title       = "Pan-Cancer — Druggable Gene Families",
    path        = file.path(out_fig, "pancancer_gene_families.png")
  )

  # n_cancers vs n_drugs scatter for druggable genes
  p_scatter <- priority_pc |>
    filter(!is.na(n_cancers)) |>
    ggplot(aes(x = n_cancers, y = n_drugs, color = Direction, size = is_hub)) +
    geom_jitter(width = 0.15, alpha = 0.75) +
    geom_text(
      data = ~ filter(., n_drugs >= quantile(n_drugs, 0.90, na.rm = TRUE) | is_hub),
      aes(label = Symbol), size = 2.8, vjust = -0.8, show.legend = FALSE
    ) +
    scale_color_manual(values = c(Up = "#de2d26", Down = "#2171b5")) +
    scale_size_manual(values = c("TRUE" = 4, "FALSE" = 2), name = "Hub gene") +
    scale_x_continuous(breaks = seq(MIN_CANCERS_PC, length(CANCERS))) +
    labs(
      title = "Pan-Cancer Druggable Targets: Recurrence vs Drug Coverage",
      x     = "Number of Cancers (significant)",
      y     = "Number of Drugs"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(out_fig, "pancancer_recurrence_vs_drugs.png"),
         plot = p_scatter, width = 10, height = 7, dpi = 300)

} else {
  message("No DGIdb results for pan-cancer DEGs.")
}

# ============================================================================
# 3. Cross-cancer comparison of druggable targets
# ============================================================================

message("Loading per-cancer drug target results for cross-cancer comparison...")

per_cancer_targets <- lapply(setNames(CANCERS, CANCERS), function(cancer) {
  path <- paste0("outputs/", cancer, "/drug_targets/", cancer, "_drug_targets.csv")
  if (!file.exists(path)) return(NULL)
  import(path) |> mutate(cancer = cancer)
})
per_cancer_targets <- Filter(Negate(is.null), per_cancer_targets)

if (length(per_cancer_targets) >= 2) {

  all_targets <- bind_rows(per_cancer_targets)

  # Genes druggable in multiple cancers
  shared_druggable <- all_targets |>
    group_by(Symbol) |>
    summarise(
      n_cancers_druggable = n(),
      cancer_list         = paste(cancer, collapse = "|"),
      mean_n_drugs        = mean(n_drugs, na.rm = TRUE),
      mean_lfc            = mean(pooled_log2FC, na.rm = TRUE),
      gene_family         = first(Gene_Family),
      any_hub             = any(is_hub),
      .groups = "drop"
    ) |>
    filter(n_cancers_druggable >= 2) |>
    arrange(desc(n_cancers_druggable), desc(mean_n_drugs))

  safe_export(shared_druggable,
              file.path(out_csv, "pancancer_shared_druggable_genes.csv"))
  message("Shared druggable genes (>=2 cancers): ", nrow(shared_druggable))

  # Drug presence matrix (gene × cancer, value = n_drugs)
  drug_matrix <- all_targets |>
    select(Symbol, cancer, n_drugs) |>
    pivot_wider(names_from = cancer, values_from = n_drugs, values_fill = 0)
  safe_export(drug_matrix, file.path(out_csv, "pancancer_drug_count_matrix.csv"))

  # ── Bubble chart: shared druggable genes × cancer
  if (nrow(shared_druggable) > 0) {
    top_shared <- head(shared_druggable, 40)
    plot_df <- all_targets |>
      filter(Symbol %in% top_shared$Symbol) |>
      mutate(Direction = ifelse(pooled_log2FC > 0, "Up", "Down"))

    p_bubble <- ggplot(plot_df,
                       aes(x = cancer, y = reorder(Symbol, n_drugs),
                           size = n_drugs, color = Direction)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c(Up = "#de2d26", Down = "#2171b5")) +
      scale_size_continuous(range = c(2, 10), name = "# Drugs") +
      labs(
        title = "Shared Druggable Targets Across Cancers",
        x     = "Cancer Type",
        y     = "Gene Symbol"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x  = element_text(angle = 35, hjust = 1),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_line(color = "gray90")
      )
    ggsave(file.path(out_fig, "pancancer_shared_druggable_bubble.png"),
           plot = p_bubble,
           width  = 10,
           height = min(48, max(7, nrow(top_shared) * 0.35 + 3)),
           dpi    = 300)
  }

  # ── Druggable gene family comparison across cancers
  family_per_cancer <- all_targets |>
    filter(Gene_Family != "Other") |>
    count(cancer, Gene_Family)

  if (nrow(family_per_cancer) > 0) {
    p_fam <- ggplot(family_per_cancer,
                    aes(x = cancer, y = n, fill = Gene_Family)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(
        title = "Druggable Gene Families per Cancer",
        x     = NULL,
        y     = "Number of Genes",
        fill  = "Gene Family"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x  = element_text(angle = 35, hjust = 1),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
      )
    ggsave(file.path(out_fig, "pancancer_family_stacked.png"),
           plot = p_fam, width = 10, height = 7, dpi = 300)
  }

  # ── Summary table: n druggable genes per cancer
  summary_tbl <- all_targets |>
    group_by(cancer) |>
    summarise(
      n_druggable   = n(),
      n_high_prio   = sum(Priority == "High"),
      n_hub_targets = sum(is_hub),
      top_target    = first(Symbol),
      .groups = "drop"
    )
  safe_export(summary_tbl, file.path(out_csv, "pancancer_drug_summary.csv"))
  message("Summary table saved.")
  print(summary_tbl)
}

message("\nPan-cancer drug target analysis complete.")
