library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)
library(msigdbr)
library(tidyverse)
library(rio)

# Create output directories
dir.create("outputs/Enrichment", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/Enrichment", showWarnings = FALSE, recursive = TRUE)

# Step 1: Load significant genes from meta-analysis
sig_genes <- import("outputs/MetaVolcanoR/significant_genes.csv")

# Step 2: Separate genes by direction
upregulated <- sig_genes |>
  filter(randomSummary > 0) |>
  pull(Gene_ID)

downregulated <- sig_genes |>
  filter(randomSummary < 0) |>
  pull(Gene_ID)

all_degs <- sig_genes |>
  pull(Gene_ID)

# Step 3: Load all genes for universe
all_genes_data <- import("outputs/MetaVolcanoR/REM.csv")
all_genes <- all_genes_data |>
  pull(Gene_ID)

# Step 4: Convert Ensembl IDs to Entrez IDs
entrez_up <- bitr(
  geneID = upregulated,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

entrez_down <- bitr(
  geneID = downregulated,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

entrez_all <- bitr(
  geneID = all_degs,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

universe <- bitr(
  geneID = all_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# ============================================================================
# GO Enrichment Analysis
# ============================================================================

# Step 5: GO enrichment for upregulated genes
go_up_bp <- enrichGO(
  gene = entrez_up$ENTREZID,
  universe = universe$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Step 6: GO enrichment for downregulated genes
go_down_bp <- enrichGO(
  gene = entrez_down$ENTREZID,
  universe = universe$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Step 7: GO enrichment for all DEGs
go_all_bp <- enrichGO(
  gene = entrez_all$ENTREZID,
  universe = universe$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# ============================================================================
# KEGG Pathway Enrichment
# ============================================================================

# Step 8: KEGG enrichment for upregulated genes
kegg_up <- enrichKEGG(
  gene = entrez_up$ENTREZID,
  universe = universe$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Step 9: KEGG enrichment for downregulated genes
kegg_down <- enrichKEGG(
  gene = entrez_down$ENTREZID,
  universe = universe$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Step 10: KEGG enrichment for all DEGs
kegg_all <- enrichKEGG(
  gene = entrez_all$ENTREZID,
  universe = universe$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# ============================================================================
# Reactome Pathway Enrichment
# ============================================================================

# Step 11: Reactome enrichment for upregulated genes
reactome_up <- enrichPathway(
  gene = entrez_up$ENTREZID,
  universe = universe$ENTREZID,
  pvalueCutoff = 0.05
)

# Step 12: Reactome enrichment for downregulated genes
reactome_down <- enrichPathway(
  gene = entrez_down$ENTREZID,
  universe = universe$ENTREZID,
  pvalueCutoff = 0.05
)

# Step 13: Reactome enrichment for all DEGs
reactome_all <- enrichPathway(
  gene = entrez_all$ENTREZID,
  universe = universe$ENTREZID,
  pvalueCutoff = 0.05
)

# ============================================================================
# Gene Set Enrichment Analysis (GSEA)
# ============================================================================

# Step 14: Prepare ranked gene list for GSEA
gene_list_data <- sig_genes |>
  left_join(entrez_all, by = c("Gene_ID" = "ENSEMBL")) |>
  filter(!is.na(ENTREZID)) |>
  select(ENTREZID, randomSummary) |>
  arrange(desc(randomSummary))

# Step 15: Create named vector for GSEA
ranked_list <- gene_list_data$randomSummary
names(ranked_list) <- gene_list_data$ENTREZID

# Step 16: GSEA for GO Biological Process
gsea_go <- gseGO(
  geneList = ranked_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05
)

# Step 17: GSEA for KEGG pathways
gsea_kegg <- gseKEGG(
  geneList = ranked_list,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05
)

# Step 18: Get MSigDB Hallmark gene sets
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") |>
  select(gs_name, entrez_gene)

# Step 19: GSEA for Hallmark gene sets
gsea_hallmark <- GSEA(
  geneList = ranked_list,
  TERM2GENE = hallmark_sets,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05
)

# ============================================================================
# Export Results
# ============================================================================

# Step 20: Export GO results
go_all_bp |>
  as.data.frame() |>
  export("outputs/Enrichment/GO_BP_results.csv")

# Step 21: Export KEGG results
kegg_all |>
  as.data.frame() |>
  export("outputs/Enrichment/KEGG_results.csv")

# Step 22: Export Reactome results
reactome_all |>
  as.data.frame() |>
  export("outputs/Enrichment/Reactome_results.csv")

# Step 23: Export GSEA GO results
gsea_go |>
  as.data.frame() |>
  export("outputs/Enrichment/GSEA_GO_results.csv")

# Step 24: Export GSEA KEGG results
gsea_kegg |>
  as.data.frame() |>
  export("outputs/Enrichment/GSEA_KEGG_results.csv")

# Step 25: Export GSEA Hallmark results
gsea_hallmark |>
  as.data.frame() |>
  export("outputs/Enrichment/GSEA_Hallmark_results.csv")

# ============================================================================
# Generate Figures
# ============================================================================

# Step 26: GO dotplot
if (nrow(as.data.frame(go_all_bp)) > 0) {
  p_go <- dotplot(go_all_bp, showCategory = 20) +
    ggtitle("GO Biological Process")

  ggsave(
    filename = "figures/Enrichment/GO_BP_dotplot.png",
    plot = p_go,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# Step 27: KEGG dotplot
if (nrow(as.data.frame(kegg_all)) > 0) {
  p_kegg <- dotplot(kegg_all, showCategory = 15) +
    ggtitle("KEGG Pathways")

  ggsave(
    filename = "figures/Enrichment/KEGG_dotplot.png",
    plot = p_kegg,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# Step 28: Reactome dotplot
if (nrow(as.data.frame(reactome_all)) > 0) {
  p_reactome <- dotplot(reactome_all, showCategory = 15) +
    ggtitle("Reactome Pathways")

  ggsave(
    filename = "figures/Enrichment/Reactome_dotplot.png",
    plot = p_reactome,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# Step 29: GSEA dotplot
if (nrow(as.data.frame(gsea_go)) > 0) {
  p_gsea <- dotplot(gsea_go, showCategory = 20) +
    ggtitle("GSEA GO Biological Process")

  ggsave(
    filename = "figures/Enrichment/GSEA_GO_dotplot.png",
    plot = p_gsea,
    width = 10,
    height = 8,
    dpi = 300
  )

  # Step 30: GSEA enrichment plot for top pathway
  top_pathway <- as.data.frame(gsea_go)$ID[1]
  p_gsea_plot <- gseaplot2(
    gsea_go,
    geneSetID = top_pathway,
    title = top_pathway
  )

  ggsave(
    filename = "figures/Enrichment/GSEA_top_pathway.png",
    plot = p_gsea_plot,
    width = 10,
    height = 6,
    dpi = 300
  )
}
