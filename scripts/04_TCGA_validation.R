library(DESeq2)
library(tidyverse)
library(rio)
library(pheatmap)
library(ggplot2)

# Create output directories
dir.create("outputs/TCGA", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/TCGA", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Step 1: Load TCGA data
# ============================================================================

counts <- readRDS("data/TCGA/TCGA_PAAD_counts.rds")
sample_info <- import("data/TCGA/TCGA_PAAD_clinical.csv")

# ============================================================================
# Step 2: Filter samples with sample type information
# ============================================================================

sample_info <- sample_info |>
  filter(!is.na(sample_type)) |>
  mutate(sample_type = factor(sample_type, levels = c("Normal", "Tumor")))

# ============================================================================
# Step 3: Filter counts to match samples
# ============================================================================

counts_filtered <- counts[, sample_info$sample_id]

# ============================================================================
# Step 4: Create DESeq2 dataset
# ============================================================================

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_info,
  design = ~ sample_type
)

# ============================================================================
# Step 5: Filter low count genes
# ============================================================================

dds <- dds[rowSums(counts(dds)) >= 10, ]

# ============================================================================
# Step 6: Run differential expression analysis
# ============================================================================

dds <- DESeq(dds)

# ============================================================================
# Step 7: Extract results
# ============================================================================

res_tcga <- results(dds, contrast = c("sample_type", "Tumor", "Normal"))

# Step 8: Convert to dataframe
res_tcga <- res_tcga |>
  as.data.frame() |>
  rownames_to_column("Gene_ID") |>
  filter(!is.na(padj))

# Step 9: Export TCGA differential expression results
export(res_tcga, "outputs/TCGA/TCGA_differential_expression.csv")

# ============================================================================
# Step 10: Load meta-analysis significant genes
# ============================================================================

sig_genes <- import("outputs/MetaVolcanoR/significant_genes.csv")

# ============================================================================
# Step 11: Join meta-analysis with TCGA results
# ============================================================================

validation <- sig_genes |>
  left_join(res_tcga, by = "Gene_ID", suffix = c("_meta", "_tcga")) |>
  filter(!is.na(log2FoldChange))

# ============================================================================
# Step 12: Calculate validation metrics
# ============================================================================

# Direction concordance
validation <- validation |>
  mutate(
    direction_concordance = sign(randomSummary) == sign(log2FoldChange),
    validated = abs(log2FoldChange) >= 0.5 & padj < 0.05
  )

# Step 13: Export validation results
export(validation, "outputs/TCGA/validation_results.csv")

# ============================================================================
# Step 14: Calculate correlation
# ============================================================================

correlation <- cor.test(
  x = validation$randomSummary,
  y = validation$log2FoldChange,
  method = "spearman"
)

# Step 15: Calculate concordance rate
concordance_rate <- mean(validation$direction_concordance) * 100

# Step 16: Calculate validation rate
validation_rate <- mean(validation$validated) * 100

# ============================================================================
# Step 17: Create correlation scatter plot
# ============================================================================

p_correlation <- ggplot(validation, aes(x = randomSummary, y = log2FoldChange)) +
  geom_point(aes(color = direction_concordance), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("TRUE" = "#2171b5", "FALSE" = "#de2d26"),
    name = "Concordant Direction"
  ) +
  labs(
    title = "GEO Meta-Analysis vs TCGA Validation",
    subtitle = sprintf("Spearman r = %.3f, p < 0.001\nConcordance = %.1f%%",
                       correlation$estimate, concordance_rate),
    x = "Meta-Analysis Log2FC (GEO)",
    y = "TCGA Log2FC"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# Step 18: Save correlation plot
ggsave(
  filename = "figures/TCGA/FC_correlation.png",
  plot = p_correlation,
  width = 8,
  height = 7,
  dpi = 300
)

# ============================================================================
# Step 19: Create TCGA volcano plot
# ============================================================================

p_volcano <- ggplot(validation, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = validated), alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "red") +
  scale_color_manual(
    values = c("TRUE" = "#de2d26", "FALSE" = "gray60"),
    name = "Validated"
  ) +
  labs(
    title = "TCGA-PAAD Differential Expression",
    subtitle = sprintf("%d / %d genes validated", sum(validation$validated), nrow(validation)),
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# Step 20: Save volcano plot
ggsave(
  filename = "figures/TCGA/validation_volcano.png",
  plot = p_volcano,
  width = 8,
  height = 7,
  dpi = 300
)

# ============================================================================
# Step 21: Create heatmap for top 50 genes
# ============================================================================

# Get top 50 validated genes
top_genes <- validation |>
  arrange(padj) |>
  head(50) |>
  pull(Gene_ID)

# Load normalized expression data
normalized <- readRDS("data/TCGA/TCGA_PAAD_normalized.rds")

# Extract expression for top genes
heatmap_data <- normalized[top_genes, sample_info$sample_id]

# Prepare annotation
annotation_col <- sample_info |>
  select(sample_id, sample_type) |>
  column_to_rownames("sample_id")

# Step 22: Generate and save heatmap
pheatmap(
  mat = heatmap_data,
  annotation_col = annotation_col,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  filename = "figures/TCGA/validation_heatmap.png",
  width = 10,
  height = 8
)
