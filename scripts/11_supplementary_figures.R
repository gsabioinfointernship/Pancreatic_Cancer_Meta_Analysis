library(tidyverse)
library(rio)
library(patchwork)
library(ggplot2)

dir.create("figures/Supplementary", showWarnings = FALSE, recursive = TRUE)

pca_dirs <- list.dirs("figures/PCA", full.names = TRUE, recursive = FALSE)

if(length(pca_dirs) > 0) {

  pca_files <- list.files("figures/PCA", pattern = "\\.png$", recursive = TRUE, full.names = TRUE)

  if(length(pca_files) > 0) {
    message("PCA plots already exist in figures/PCA/")
    message("Total PCA plots: ", length(pca_files))
  }
}

geo_ids <- c("GSE130688", "GSE136569", "GSE171485", "GSE196009",
             "GSE211398", "GSE280271", "GSE293744")

for(geo_id in geo_ids) {
  if(file.exists(paste0("outputs/DESeq2/", geo_id, ".csv"))) {

    deseq_results <- import(paste0("outputs/DESeq2/", geo_id, ".csv"))

    p_volcano <- ggplot(deseq_results, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = abs(log2FoldChange) >= 1 & padj < 0.05), alpha = 0.6, size = 1.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
      scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#de2d26"),
                         name = "Significant") +
      labs(
        title = paste0(geo_id, " Differential Expression"),
        x = "Log2 Fold Change",
        y = "-Log10(Adjusted P-value)"
      ) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    ggsave(paste0("figures/Supplementary/", geo_id, "_volcano.png"),
           p_volcano, width = 8, height = 6, dpi = 300)
  }
}

if(file.exists("outputs/Enrichment/Reactome_results.csv")) {

  reactome <- import("outputs/Enrichment/Reactome_results.csv") %>%
    arrange(p.adjust) %>%
    head(15)

  if(nrow(reactome) > 0) {
    p_reactome <- ggplot(reactome, aes(x = Count, y = reorder(Description, -p.adjust))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "#de2d26", high = "#2171b5", name = "Adj. P-value") +
      scale_size_continuous(name = "Gene Count") +
      labs(title = "Reactome Pathways", x = "Gene Count", y = "") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    ggsave("figures/Supplementary/SuppFig_Reactome.png", p_reactome,
           width = 10, height = 6, dpi = 300)
  }
}

if(file.exists("outputs/Survival/Cox_univariate.csv")) {

  cox_results <- import("outputs/Survival/Cox_univariate.csv") %>%
    filter(FDR < 0.05) %>%
    arrange(p_value) %>%
    head(30)

  if(nrow(cox_results) > 0) {
    cox_results <- cox_results %>%
      mutate(Gene_ID = factor(Gene_ID, levels = rev(Gene_ID)))

    p_forest <- ggplot(cox_results, aes(x = HR, y = Gene_ID)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.3) +
      geom_point(aes(color = FDR < 0.01), size = 3) +
      scale_color_manual(values = c("TRUE" = "#de2d26", "FALSE" = "#2171b5")) +
      scale_x_continuous(trans = "log10") +
      labs(
        title = "Univariate Cox Regression (Top 30 Genes)",
        x = "Hazard Ratio (95% CI, log scale)",
        y = "Gene"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none"
      )

    ggsave("figures/Supplementary/SuppFig_Forest_Univariate.png", p_forest,
           width = 10, height = 10, dpi = 300)
  }
}

if(file.exists("outputs/MachineLearning/model_performance.csv")) {

  ml_all <- import("outputs/MachineLearning/model_performance.csv")

  ml_validation <- ml_all %>%
    filter(Dataset == "Validation (TCGA)")

  p_ml_comparison <- ggplot(ml_validation, aes(x = Model, y = AUC, fill = Feature_Set)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = round(AUC, 3)), position = position_dodge(width = 0.9),
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("top50" = "#2171b5", "top100" = "#fdae61", "top200" = "#de2d26"),
                      name = "Feature Set") +
    ylim(0, 1) +
    labs(
      title = "ML Performance Across Feature Sets",
      x = "Model",
      y = "AUC-ROC (TCGA Validation)"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("figures/Supplementary/SuppFig_ML_FeatureSets.png", p_ml_comparison,
         width = 10, height = 6, dpi = 300)

  ml_metrics <- ml_validation %>%
    filter(Feature_Set == "top100") %>%
    pivot_longer(cols = c(AUC, Accuracy, Sensitivity, Specificity),
                 names_to = "Metric", values_to = "Value")

  p_metrics <- ggplot(ml_metrics, aes(x = Model, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = "Set2") +
    ylim(0, 1) +
    labs(
      title = "ML Performance Metrics (Top 100 Features)",
      x = "Model",
      y = "Value"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("figures/Supplementary/SuppFig_ML_Metrics.png", p_metrics,
         width = 10, height = 6, dpi = 300)
}

if(file.exists("outputs/Network/network_modules.csv")) {

  modules <- import("outputs/Network/network_modules.csv")

  module_summary <- modules %>%
    group_by(Module) %>%
    summarise(
      n_genes = n(),
      avg_degree = mean(Degree),
      n_upregulated = sum(Log2FC > 0),
      n_downregulated = sum(Log2FC < 0)
    ) %>%
    arrange(desc(n_genes)) %>%
    head(10)

  p_modules <- ggplot(module_summary, aes(x = reorder(as.factor(Module), n_genes), y = n_genes)) +
    geom_bar(stat = "identity", fill = "#2171b5") +
    coord_flip() +
    labs(
      title = "Top 10 Network Modules by Size",
      x = "Module ID",
      y = "Number of Genes"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("figures/Supplementary/SuppFig_Network_Modules.png", p_modules,
         width = 8, height = 6, dpi = 300)
}

if(file.exists("outputs/DrugTarget/druggable_genes.csv")) {

  druggable <- import("outputs/DrugTarget/druggable_genes.csv")

  if("druggable_class" %in% colnames(druggable)) {
    class_count <- druggable %>%
      group_by(druggable_class) %>%
      summarise(count = n()) %>%
      arrange(desc(count))

    p_drug_class <- ggplot(class_count, aes(x = reorder(druggable_class, count), y = count)) +
      geom_bar(stat = "identity", fill = "#de2d26") +
      coord_flip() +
      labs(
        title = "Druggable Gene Classes",
        x = "Drug Target Class",
        y = "Number of Genes"
      ) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    ggsave("figures/Supplementary/SuppFig_DrugClasses.png", p_drug_class,
           width = 8, height = 6, dpi = 300)
  }
}

message("Supplementary figures created successfully!")
message("Check figures/Supplementary/ for all outputs")
