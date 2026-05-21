source("drug_target_analysis/00_drug_functions.R")

for (cancer in CANCERS) {

  message("\n=== ", toupper(cancer), " ===")

  # ── Directories
  out_csv <- paste0("outputs/", cancer, "/drug_targets/")
  out_fig <- paste0("outputs/", cancer, "/figures/drug_targets/")
  dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

  # ── Load data
  meta <- tryCatch(load_meta(cancer), error = function(e) {
    message("  SKIPPING — cannot load meta: ", e$message); NULL
  })
  if (is.null(meta)) next

  hub_df <- load_hub_genes(cancer)
  if (is.null(hub_df)) {
    message("  Hub genes not found — all sig DEGs used as candidates")
  } else {
    message("  Hub genes loaded: ", nrow(hub_df))
  }

  # ── Significant DEGs as query genes
  sig <- meta |>
    filter(rem_padj < PADJ_CUTOFF, abs(pooled_log2FC) >= LFC_CUTOFF) |>
    filter(!is.na(Symbol), Symbol != "")

  message("  Significant DEGs: ", nrow(sig))

  if (nrow(sig) == 0) {
    message("  SKIPPING — no significant DEGs")
    next
  }

  # ── DGIdb query
  message("  Querying DGIdb...")
  iact <- query_dgidb(sig$Symbol)

  if (is.null(iact)) {
    message("  No DGIdb results — exporting sig DEG list only")
    safe_export(sig, file.path(out_csv, paste0(cancer, "_significant_DEGs.csv")))
    next
  }

  message("  Drug interactions found: ", nrow(iact))
  message("  Druggable genes: ", n_distinct(iact$Symbol))

  # ── Priority table
  priority <- build_priority_table(iact, meta, hub_df)

  # ── Export
  export_drug(iact, priority, out_csv, cancer)

  # ── Summary stats
  high_priority <- filter(priority, Priority == "High")
  message("  High-priority targets (hub + >=3 drugs): ", nrow(high_priority))

  if (nrow(high_priority) > 0) {
    safe_export(high_priority,
                file.path(out_csv, paste0(cancer, "_high_priority_targets.csv")))
  }

  # ── Figures
  message("  Saving figures...")

  save_drug_lollipop(
    priority_df = priority,
    title       = paste(tools::toTitleCase(cancer), "— Top Druggable Targets"),
    path        = file.path(out_fig, paste0(cancer, "_drug_lollipop.png"))
  )

  save_family_barplot(
    priority_df = priority,
    title       = paste(tools::toTitleCase(cancer), "— Druggable Gene Families"),
    path        = file.path(out_fig, paste0(cancer, "_gene_families.png"))
  )

  if (!is.null(hub_df)) {
    save_hub_drug_boxplot(
      priority_df = priority,
      title       = paste(tools::toTitleCase(cancer), "— Hub vs Non-Hub Drug Count"),
      path        = file.path(out_fig, paste0(cancer, "_hub_drug_boxplot.png"))
    )
  }

  # Direction breakdown of druggable genes
  dir_summary <- priority |>
    count(Direction, Gene_Family) |>
    filter(Gene_Family != "Other")

  if (nrow(dir_summary) > 0) {
    p_dir <- ggplot(dir_summary,
                    aes(x = reorder(Gene_Family, n), y = n, fill = Direction)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c(Up = "#de2d26", Down = "#2171b5")) +
      coord_flip() +
      labs(
        title = paste(tools::toTitleCase(cancer), "— Druggable Families by Direction"),
        x = NULL, y = "Count"
      ) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave(file.path(out_fig, paste0(cancer, "_families_direction.png")),
           plot = p_dir, width = 9, height = 6, dpi = 300)
  }

  message("  Done: ", cancer)
}

message("\nPer-cancer drug target analysis complete.")
