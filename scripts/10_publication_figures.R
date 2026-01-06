library(tidyverse)
library(rio)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)

dir.create("figures/Publication", showWarnings = FALSE, recursive = TRUE)

sig_genes <- import("outputs/MetaVolcanoR/significant_genes.csv")

p_volcano <- ggplot(sig_genes, aes(x = randomSummary, y = -log10(randomP))) +
  geom_point(aes(color = abs(randomSummary) >= 1 & randomP < 0.05), alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#de2d26"),
                     name = "Significant") +
  labs(
    title = "Meta-Analysis: 613 Significant DEGs",
    x = "Log2 Fold Change (Random Effects)",
    y = "-Log10(P-value)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "right"
  )

top_50_genes <- sig_genes %>%
  arrange(randomP) %>%
  head(50) %>%
  pull(Gene_ID)

if(file.exists("data/TCGA/TCGA_PAAD_normalized.rds") &&
   file.exists("data/TCGA/TCGA_PAAD_clinical.csv")) {

  normalized <- readRDS("data/TCGA/TCGA_PAAD_normalized.rds")
  clinical <- import("data/TCGA/TCGA_PAAD_clinical.csv")

  heatmap_data <- normalized[top_50_genes, clinical$sample_id]

  annotation_col <- clinical %>%
    dplyr::select(sample_id, sample_type) %>%
    column_to_rownames("sample_id") %>%
    as.data.frame()

  col_fun <- colorRamp2(c(-3, 0, 3), c("#2171b5", "white", "#de2d26"))

  ht <- Heatmap(
    heatmap_data,
    name = "Expression",
    col = col_fun,
    top_annotation = HeatmapAnnotation(
      Type = annotation_col$sample_type,
      col = list(Type = c("Normal" = "#2171b5", "Tumor" = "#de2d26"))
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    column_title = "Top 50 DEGs in TCGA-PAAD"
  )

  png("figures/Publication/Figure2_heatmap.png", width = 10, height = 8, units = "in", res = 300)
  draw(ht)
  dev.off()
}

ggsave("figures/Publication/Figure2_volcano.png", p_volcano, width = 10, height = 8, dpi = 300)

if(file.exists("outputs/Enrichment/GO_BP_results.csv") &&
   file.exists("outputs/Enrichment/KEGG_results.csv")) {

  go_results <- import("outputs/Enrichment/GO_BP_results.csv") %>%
    arrange(p.adjust) %>%
    head(15)

  kegg_results <- import("outputs/Enrichment/KEGG_results.csv") %>%
    arrange(p.adjust) %>%
    head(10)

  if(nrow(go_results) > 0) {
    p_go <- ggplot(go_results, aes(x = Count, y = reorder(Description, -p.adjust))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "#de2d26", high = "#2171b5", name = "Adj. P-value") +
      scale_size_continuous(name = "Gene Count") +
      labs(title = "GO Biological Process", x = "Gene Count", y = "") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    ggsave("figures/Publication/Figure3_GO.png", p_go, width = 10, height = 6, dpi = 300)
  }

  if(nrow(kegg_results) > 0) {
    p_kegg <- ggplot(kegg_results, aes(x = Count, y = reorder(Description, -p.adjust))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "#de2d26", high = "#2171b5", name = "Adj. P-value") +
      scale_size_continuous(name = "Gene Count") +
      labs(title = "KEGG Pathways", x = "Gene Count", y = "") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    ggsave("figures/Publication/Figure3_KEGG.png", p_kegg, width = 10, height = 6, dpi = 300)
  }
}

if(file.exists("figures/TCGA/FC_correlation.png") &&
   file.exists("figures/TCGA/validation_heatmap.png")) {

  message("TCGA validation figures already created in script 04")
}

if(file.exists("outputs/MachineLearning/model_performance.csv")) {

  ml_results <- import("outputs/MachineLearning/model_performance.csv") %>%
    filter(Dataset == "Validation (TCGA)", Feature_Set == "top100")

  p_ml <- ggplot(ml_results, aes(x = Model, y = AUC)) +
    geom_bar(stat = "identity", fill = "#2171b5", width = 0.7) +
    geom_text(aes(label = round(AUC, 3)), vjust = -0.5, size = 4) +
    ylim(0, 1) +
    labs(
      title = "Machine Learning Performance (TCGA Validation)",
      x = "Model",
      y = "AUC-ROC"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("figures/Publication/Figure6_ML_performance.png", p_ml, width = 8, height = 6, dpi = 300)
}

dataset_summary <- tibble(
  Dataset = c("GSE130688", "GSE136569", "GSE171485", "GSE196009",
              "GSE211398", "GSE280271", "GSE293744", "Total"),
  Normal = c(15, 5, 6, 6, 12, 7, 8, 59),
  Tumor = c(15, 5, 6, 13, 16, 6, 40, 101),
  Total = c(30, 10, 12, 19, 28, 13, 48, 160)
)

export(dataset_summary, "figures/Publication/Table1_dataset_summary.csv")

message("Main publication figures created successfully!")
message("Check figures/Publication/ for all outputs")
