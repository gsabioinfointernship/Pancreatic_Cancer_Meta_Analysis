source("regulatory_analysis/00_regulatory_functions.R")
library(pheatmap)

# 1. Setup
pancancer_csv <- "outputs/pancancer/regulatory_analysis/"
pancancer_fig <- "outputs/pancancer/figures/regulatory_analysis/"
dir.create(pancancer_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(pancancer_fig, showWarnings = FALSE, recursive = TRUE)

message("\n=== PAN-CANCER REGULATORY ANALYSIS AGGREGATION ===")

# ── 1. Aggregate TF Enrichment (UP)
up_files <- list.files("outputs", pattern = "_TF_enrichment_up.csv", recursive = TRUE, full.names = TRUE)

if (length(up_files) > 0) {
  message("  Aggregating UP-regulated TF enrichment...")
  all_tf_up <- lapply(up_files, function(f) {
    cancer_name <- gsub("_TF_enrichment_up.csv", "", basename(f))
    read.csv(f) |> 
      mutate(
        cancer = cancer_name,
        across(any_of(c("ID", "Description", "GeneRatio", "BgRatio", "geneID")), as.character)
      )
  }) |> bind_rows()
  
  write.csv(all_tf_up, file.path(pancancer_csv, "pancancer_all_TF_up.csv"), row.names = FALSE)
  
  # Identify universal master regulators
  summary_tf_up <- all_tf_up |>
    group_by(ID) |>
    summarise(
      n_cancers = n_distinct(cancer),
      cancers = paste(unique(cancer), collapse = "; "),
      mean_qval = mean(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(n_cancers), mean_qval)
    
  write.csv(summary_tf_up, file.path(pancancer_csv, "pancancer_consistent_TFs_up.csv"), row.names = FALSE)
  
  # Heatmap of consistent TFs
  top_tfs_up <- summary_tf_up |> filter(n_cancers >= 2) |> head(30) |> pull(ID)
  
  if (length(top_tfs_up) >= 2) {
    plot_df <- all_tf_up |> 
      filter(ID %in% top_tfs_up) |>
      select(ID, cancer, p.adjust) |>
      mutate(logP = -log10(p.adjust)) |>
      pivot_wider(id_cols = ID, names_from = cancer, values_from = logP, values_fill = 0) |>
      column_to_rownames("ID")
      
    pheatmap(
      as.matrix(plot_df),
      color = colorRampPalette(c("white", "#de2d26"))(100),
      main = "Consistent Master Regulators (UP, -log10 Adj.P)",
      filename = file.path(pancancer_fig, "pancancer_consistent_TF_up_heatmap.png"),
      width = 8, height = max(6, length(top_tfs_up) * 0.3)
    )
  }
}

# ── 2. Aggregate TF Enrichment (DOWN)
down_files <- list.files("outputs", pattern = "_TF_enrichment_down.csv", recursive = TRUE, full.names = TRUE)

if (length(down_files) > 0) {
  message("  Aggregating DOWN-regulated TF enrichment...")
  all_tf_down <- lapply(down_files, function(f) {
    cancer_name <- gsub("_TF_enrichment_down.csv", "", basename(f))
    read.csv(f) |> 
      mutate(
        cancer = cancer_name,
        across(any_of(c("ID", "Description", "GeneRatio", "BgRatio", "geneID")), as.character)
      )
  }) |> bind_rows()
  
  write.csv(all_tf_down, file.path(pancancer_csv, "pancancer_all_TF_down.csv"), row.names = FALSE)
  
  summary_tf_down <- all_tf_down |>
    group_by(ID) |>
    summarise(
      n_cancers = n_distinct(cancer),
      cancers = paste(unique(cancer), collapse = "; "),
      mean_qval = mean(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(n_cancers), mean_qval)
    
  write.csv(summary_tf_down, file.path(pancancer_csv, "pancancer_consistent_TFs_down.csv"), row.names = FALSE)
  
  top_tfs_down <- summary_tf_down |> filter(n_cancers >= 2) |> head(30) |> pull(ID)
  
  if (length(top_tfs_down) >= 2) {
    plot_df <- all_tf_down |> 
      filter(ID %in% top_tfs_down) |>
      select(ID, cancer, p.adjust) |>
      mutate(logP = -log10(p.adjust)) |>
      pivot_wider(id_cols = ID, names_from = cancer, values_from = logP, values_fill = 0) |>
      column_to_rownames("ID")
      
    pheatmap(
      as.matrix(plot_df),
      color = colorRampPalette(c("white", "#3182bd"))(100),
      main = "Consistent Master Regulators (DOWN, -log10 Adj.P)",
      filename = file.path(pancancer_fig, "pancancer_consistent_TF_down_heatmap.png"),
      width = 8, height = max(6, length(top_tfs_down) * 0.3)
    )
  }
}

# ── 3. Aggregate TF GSEA
gsea_files <- list.files("outputs", pattern = "_TF_GSEA.csv", recursive = TRUE, full.names = TRUE)

if (length(gsea_files) > 0) {
  message("  Aggregating TF GSEA results...")
  all_tf_gsea <- lapply(gsea_files, function(f) {
    cancer_name <- gsub("_TF_GSEA.csv", "", basename(f))
    read.csv(f) |> 
      mutate(
        cancer = cancer_name,
        across(any_of(c("ID", "Description", "core_enrichment")), as.character)
      )
  }) |> bind_rows()
  
  write.csv(all_tf_gsea, file.path(pancancer_csv, "pancancer_all_TF_GSEA.csv"), row.names = FALSE)
  
  summary_gsea <- all_tf_gsea |>
    group_by(ID) |>
    summarise(
      n_cancers = n_distinct(cancer),
      cancers = paste(unique(cancer), collapse = "; "),
      mean_NES = mean(NES, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(n_cancers), desc(abs(mean_NES)))
    
  write.csv(summary_gsea, file.path(pancancer_csv, "pancancer_consistent_GSEA_TFs.csv"), row.names = FALSE)
}

message("  Pan-cancer regulatory aggregation complete.")
