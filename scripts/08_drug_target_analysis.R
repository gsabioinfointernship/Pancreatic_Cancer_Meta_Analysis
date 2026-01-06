library(tidyverse)
library(rio)

dir.create("outputs/DrugTarget", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/DrugTarget", showWarnings = FALSE, recursive = TRUE)

sig_genes <- import("outputs/MetaVolcanoR/significant_genes.csv")
hub_genes <- import("outputs/Network/hub_genes.csv")

priority_genes <- sig_genes %>%
  filter(!is.na(Gene_Symbol), Gene_Symbol != "") %>%
  mutate(is_hub = Gene_Symbol %in% hub_genes$Gene) %>%
  arrange(randomP)

if(requireNamespace("rDGIdb", quietly = TRUE)) {

  library(rDGIdb)

  gene_symbols <- priority_genes$Gene_Symbol

  drug_interactions <- queryDGIdb(
    genes = gene_symbols,
    sourceDatabases = c("DrugBank", "ChEMBL", "PharmGKB"),
    interactionTypes = NULL
  )

  if(!is.null(drug_interactions) && nrow(drug_interactions) > 0) {
    drug_data <- drug_interactions %>%
      as.data.frame() %>%
      left_join(priority_genes %>% dplyr::select(Gene_Symbol, randomSummary, randomP, is_hub),
                by = c("geneName" = "Gene_Symbol"))

    export(drug_data, "outputs/DrugTarget/drug_interactions.csv")

    druggable_genes <- drug_data %>%
      group_by(geneName) %>%
      summarise(
        n_drugs = n_distinct(drugName),
        n_interactions = n(),
        drug_names = paste(unique(drugName), collapse = "; "),
        interaction_types = paste(unique(interactionType), collapse = "; "),
        Log2FC = first(randomSummary),
        P_value = first(randomP),
        is_hub = first(is_hub)
      ) %>%
      arrange(desc(n_drugs))

    export(druggable_genes, "outputs/DrugTarget/druggable_genes.csv")

    message("Drug-gene interaction analysis completed using rDGIdb!")
    message("Druggable genes found: ", nrow(druggable_genes))
    message("Total drug interactions: ", nrow(drug_data))

  } else {
    message("No drug interactions found via rDGIdb. Creating placeholder files.")

    druggable_genes <- tibble(
      geneName = character(),
      n_drugs = numeric(),
      drug_names = character(),
      Log2FC = numeric(),
      is_hub = logical()
    )

    export(druggable_genes, "outputs/DrugTarget/druggable_genes.csv")
  }

} else {

  message("rDGIdb package not available. Creating manual druggable gene classification.")

  druggable_families <- c(
    "kinase", "receptor", "channel", "transporter", "protease",
    "phosphatase", "transferase", "hydrolase", "oxidoreductase",
    "GPCR", "nuclear receptor", "transcription factor"
  )

  druggable_genes <- priority_genes %>%
    filter(str_detect(Gene_Description, regex(paste(druggable_families, collapse = "|"), ignore_case = TRUE))) %>%
    mutate(
      druggable_class = case_when(
        str_detect(Gene_Description, regex("kinase", ignore_case = TRUE)) ~ "Kinase",
        str_detect(Gene_Description, regex("receptor", ignore_case = TRUE)) ~ "Receptor",
        str_detect(Gene_Description, regex("channel", ignore_case = TRUE)) ~ "Ion Channel",
        str_detect(Gene_Description, regex("transporter", ignore_case = TRUE)) ~ "Transporter",
        str_detect(Gene_Description, regex("protease", ignore_case = TRUE)) ~ "Protease",
        TRUE ~ "Other Druggable"
      )
    ) %>%
    dplyr::select(Gene_Symbol, Gene_Description, druggable_class, randomSummary, randomP, is_hub)

  export(druggable_genes, "outputs/DrugTarget/druggable_genes.csv")

  message("Manual druggable gene classification completed!")
  message("Druggable genes found: ", nrow(druggable_genes))
}

if(file.exists("outputs/DrugTarget/druggable_genes.csv")) {

  druggable_data <- import("outputs/DrugTarget/druggable_genes.csv")

  if(nrow(druggable_data) > 0) {

    if("druggable_class" %in% colnames(druggable_data)) {
      class_summary <- druggable_data %>%
        group_by(druggable_class) %>%
        summarise(count = n()) %>%
        arrange(desc(count))

      p1 <- ggplot(class_summary, aes(x = reorder(druggable_class, count), y = count)) +
        geom_bar(stat = "identity", fill = "#2171b5") +
        coord_flip() +
        labs(title = "Druggable Gene Classes",
             x = "Gene Class", y = "Number of Genes") +
        theme_bw() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))

      ggsave("figures/DrugTarget/druggable_categories.png", p1,
             width = 8, height = 6, dpi = 300)
    }

    if("n_drugs" %in% colnames(druggable_data)) {
      top_targets <- druggable_data %>%
        arrange(desc(n_drugs)) %>%
        head(20)

      p2 <- ggplot(top_targets, aes(x = reorder(geneName, n_drugs), y = n_drugs)) +
        geom_bar(stat = "identity", aes(fill = Log2FC)) +
        scale_fill_gradient2(low = "#2171b5", mid = "white", high = "#de2d26",
                             midpoint = 0, name = "Log2FC") +
        coord_flip() +
        labs(title = "Top 20 Druggable Targets by Number of Drugs",
             x = "Gene", y = "Number of Drugs") +
        theme_bw() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))

      ggsave("figures/DrugTarget/top_druggable_targets.png", p2,
             width = 8, height = 6, dpi = 300)
    }
  }
}

message("Drug target analysis completed!")
