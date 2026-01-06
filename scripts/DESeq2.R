# Differential Gene Expression Analysis
# Author: Muntasim Fuad

# 00. Load required packages
library(tidyverse)
library(DESeq2)
library(sva)
library(rio)


# 01. Load count data
GSE130688 <- import("data/GSE202853_raw_counts.tsv")
glimpse(GSE130688)

# 02. Load metadata
metadata <- import("data/GSE130688_metadata.csv")
glimpse(metadata)

# 03. Create a matrix & add gene ids as row names
count_data <- GSE130688 |> column_to_rownames("GeneID") |>
  as.matrix()

# 04. Match metadata with count data
metadata <- metadata |> 
  filter(sample %in% colnames(count_data)) |> 
  arrange(match(sample, colnames(count_data)))

# 05. Prepare Sample information
colData <- data.frame( condition = as.factor(metadata$condition), 
                       row.names = colnames(count_data))

# 06. Create DESeq2 data set object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = colData,
                              design = ~ condition)

# 07. filter any counts less than 10
keep <- rowSums(counts(dds) > round(nrow(metadata)/2)) >= 10
dds <- dds[keep,]

# 06. Estimate number of surrogate variables
mod    <- model.matrix(~ condition, data = colData(dds))
mod0   <- model.matrix(~ 1,         data = colData(dds))
svobj  <- svaseq(counts(dds), mod, mod0)
n.sv <- svobj$n.sv  

# 07. Add SVs to colData and update design
for(i in seq_len(n.sv)) {
  colData(dds)[[paste0("SV", i)]] <- svobj$sv[, i]
}
sv_terms <- paste0("SV", seq_len(n.sv), collapse = " + ")
design(dds) <- as.formula(paste("~", sv_terms, "+ condition"))

# 08. Run DESeq once with full design
dds_sva <- DESeq(dds)
# ------------------------------------------------------------------------------
# 09. Analyze the data set and compile result
res <- results(dds_sva)
res$Gene_ID <- row.names(res)

# 12. Export results to a CSV file
dir.create("outputs/DESeq2/", showWarnings = FALSE)
write.csv(res, "outputs/DESeq2/GSE130688.csv",
          row.names = FALSE)
# ------------------------------------------------------------------------------
# 09. Variance-stabilizing transformation for PCA
vsd_pre  <- vst(dds, blind = TRUE)
vsd_post <- vst(dds_sva, blind = FALSE)

# PCA plots (with corrected box.padding and shape)
pca_data_pre <- plotPCA(vsd_pre,  intgroup = "condition", returnData = TRUE)
pca_data_post <- plotPCA(vsd_post, intgroup =  "condition", returnData = TRUE)

percent_pre  <- round(100 * attr(pca_data_pre,  "percentVar"))
percent_post <- round(100 * attr(pca_data_post, "percentVar"))

# Common color palette
custom_colors <- c("Cancer" = "#de2d26",
                   "Normal" = "#2171b5")

# Pre-correction plot
pca_pre <- ggplot(pca_data_pre, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_pre[1], "% variance")) +
  ylab(paste0("PC2: ", percent_pre[2], "% variance")) +
  ggtitle("GSE130688 (Pre-Correction)") +
  scale_color_manual(values = custom_colors) +
  theme_minimal() + theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    panel.border  = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    aspect.ratio   = 1
  )

# Post-correction plot
pca_post <- ggplot(pca_data_post, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_post[1], "% variance")) +
  ylab(paste0("PC2: ", percent_post[2], "% variance")) +
  ggtitle("GSE130688 (Post-Correction)") +
  scale_color_manual(values = custom_colors) +
  theme_minimal() + theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    panel.border  = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    aspect.ratio   = 1
  )

# Save the PCA plots
dir.create("figures/PCA", showWarnings = FALSE)
dir.create("figures/PCA/GSE130688", showWarnings = FALSE)
ggsave(filename = "figures/PCA/GSE130688/GSE130688_pre.png",
       plot = pca_pre,
       width = 5,
       height = 5,
       units = "in",
       bg = "white",
       dpi = 600)


ggsave(filename = "figures/PCA/GSE130688/GSE130688_post.png",
       plot = pca_post,
       width = 5,
       height = 5,
       units = "in",
       bg = "white",
       dpi = 600)