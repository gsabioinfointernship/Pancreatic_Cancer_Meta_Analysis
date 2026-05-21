source("regulatory_analysis/00_regulatory_functions.R")
library(msigdbr)

# 1. Load MSigDB TF Targets (C3:TFT)
message("Loading MSigDB TF-target database...")
tf_db <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") |>
  select(gs_name, gene_symbol) |>
  rename(TF = gs_name, Target = gene_symbol)

# Clean TF names (e.g., "TF_TARGET_GENES" -> "TF")
tf_db$TF <- gsub("_TARGET_GENES", "", tf_db$TF)

for (cancer in CANCERS) {
  
  message("\n=== ", toupper(cancer), " ===")
  
  # 1. Setup
  out_csv <- file.path("outputs", cancer, "regulatory_analysis")
  out_fig <- file.path("outputs", cancer, "figures", "regulatory_analysis")
  dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
  
  # 2. Load Meta results
  meta <- load_meta(cancer)
  if (is.null(meta)) next
  
  # 3. Master Regulator ORA (Up-regulated)
  message("    Running Master Regulator Analysis (Up)...")
  tf_up <- run_tf_ora(meta, tf_db, direction = "up")
  if (!is.null(tf_up)) {
    write.csv(tf_up, file.path(out_csv, paste0(cancer, "_TF_enrichment_up.csv")), row.names = FALSE)
    save_tf_dotplot(tf_up, cancer, "Up", file.path(out_fig, paste0(cancer, "_TF_dotplot_up.png")))
  }
  
  # 4. Master Regulator ORA (Down-regulated)
  message("    Running Master Regulator Analysis (Down)...")
  tf_down <- run_tf_ora(meta, tf_db, direction = "down")
  if (!is.null(tf_down)) {
    write.csv(tf_down, file.path(out_csv, paste0(cancer, "_TF_enrichment_down.csv")), row.names = FALSE)
    save_tf_dotplot(tf_down, cancer, "Down", file.path(out_fig, paste0(cancer, "_TF_dotplot_down.png")))
  }
  
  # 5. Master Regulator GSEA (Global)
  message("    Running Master Regulator GSEA...")
  tf_gsea <- run_tf_gsea(meta, tf_db)
  if (!is.null(tf_gsea)) {
    write.csv(tf_gsea, file.path(out_csv, paste0(cancer, "_TF_GSEA.csv")), row.names = FALSE)
  }
  
  message("    Done: ", cancer)
}

message("\nPer-cancer regulatory analysis complete.")
