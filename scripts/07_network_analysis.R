library(STRINGdb)
library(igraph)
library(tidyverse)
library(rio)
library(ggraph)

dir.create("outputs/Network", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/Network", showWarnings = FALSE, recursive = TRUE)

sig_genes <- import("outputs/MetaVolcanoR/significant_genes.csv")

gene_symbols <- sig_genes %>%
  dplyr::select(Gene_ID, Gene_Symbol, randomSummary) %>%
  filter(!is.na(Gene_Symbol), Gene_Symbol != "")

string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 700)

mapped <- string_db$map(gene_symbols, "Gene_Symbol", removeUnmappedRows = TRUE)

interactions <- string_db$get_interactions(mapped$STRING_id)

ppi_network <- interactions %>%
  left_join(mapped %>% dplyr::select(STRING_id, Gene_Symbol, randomSummary),
            by = c("from" = "STRING_id")) %>%
  dplyr::rename(Gene_from = Gene_Symbol, FC_from = randomSummary) %>%
  left_join(mapped %>% dplyr::select(STRING_id, Gene_Symbol, randomSummary),
            by = c("to" = "STRING_id")) %>%
  dplyr::rename(Gene_to = Gene_Symbol, FC_to = randomSummary) %>%
  filter(!is.na(Gene_from), !is.na(Gene_to))

export(ppi_network, "outputs/Network/PPI_network.csv")

g <- graph_from_data_frame(
  d = ppi_network %>% dplyr::select(Gene_from, Gene_to, combined_score),
  vertices = mapped %>% dplyr::select(Gene_Symbol, randomSummary),
  directed = FALSE
)

g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

degree_cent <- degree(g)
betweenness_cent <- betweenness(g)
closeness_cent <- closeness(g)
eigen_cent <- eigen_centrality(g)$vector

network_metrics <- tibble(
  Gene = V(g)$name,
  Degree = degree_cent,
  Betweenness = betweenness_cent,
  Closeness = closeness_cent,
  Eigenvector = eigen_cent,
  Log2FC = V(g)$randomSummary
) %>%
  arrange(desc(Degree))

hub_threshold <- quantile(network_metrics$Degree, 0.90)
hub_genes <- network_metrics %>%
  filter(Degree >= hub_threshold) %>%
  arrange(desc(Degree))

export(hub_genes, "outputs/Network/hub_genes.csv")

communities <- cluster_louvain(g)
V(g)$community <- membership(communities)

network_modules <- tibble(
  Gene = V(g)$name,
  Module = V(g)$community,
  Degree = degree_cent,
  Log2FC = V(g)$randomSummary
) %>%
  arrange(Module, desc(Degree))

export(network_modules, "outputs/Network/network_modules.csv")

module_sizes <- table(V(g)$community)
top_modules <- names(sort(module_sizes, decreasing = TRUE)[1:5])

V(g)$size <- scales::rescale(degree_cent, to = c(3, 15))
V(g)$color <- ifelse(V(g)$randomSummary > 0, "#de2d26", "#2171b5")

png("figures/Network/PPI_network.png", width = 12, height = 12, units = "in", res = 300)
set.seed(123)
plot(g,
     vertex.size = V(g)$size,
     vertex.color = V(g)$color,
     vertex.label = ifelse(degree_cent >= hub_threshold, V(g)$name, NA),
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     edge.width = 0.5,
     edge.color = "gray70",
     layout = layout_with_fr(g),
     main = "Protein-Protein Interaction Network")
legend("topright",
       legend = c("Upregulated", "Downregulated", "Hub gene"),
       col = c("#de2d26", "#2171b5", "black"),
       pch = c(19, 19, NA),
       lty = c(NA, NA, NA),
       bty = "n")
dev.off()

p_hubs <- ggplot(hub_genes %>% head(20), aes(x = reorder(Gene, Degree), y = Degree)) +
  geom_bar(stat = "identity", aes(fill = Log2FC)) +
  scale_fill_gradient2(low = "#2171b5", mid = "white", high = "#de2d26",
                       midpoint = 0, name = "Log2FC") +
  coord_flip() +
  labs(title = "Top 20 Hub Genes by Degree Centrality",
       x = "Gene", y = "Degree Centrality") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/Network/hub_genes_barplot.png", p_hubs, width = 8, height = 6, dpi = 300)

for(module_id in top_modules) {
  module_genes <- network_modules %>% filter(Module == module_id) %>% pull(Gene)

  subg <- induced_subgraph(g, vids = which(V(g)$name %in% module_genes))

  if(vcount(subg) >= 5) {
    png(paste0("figures/Network/module_", module_id, "_network.png"),
        width = 8, height = 8, units = "in", res = 300)
    set.seed(123)
    plot(subg,
         vertex.size = scales::rescale(degree(subg), to = c(5, 15)),
         vertex.color = V(subg)$color,
         vertex.label.cex = 0.6,
         vertex.label.color = "black",
         edge.width = 1,
         edge.color = "gray70",
         layout = layout_with_fr(subg),
         main = paste0("Network Module ", module_id, " (", vcount(subg), " genes)"))
    dev.off()
  }
}

message("Network analysis completed!")
message("Total genes in network: ", vcount(g))
message("Total interactions: ", ecount(g))
message("Hub genes (top 10%): ", nrow(hub_genes))
message("Network modules detected: ", length(unique(V(g)$community)))
