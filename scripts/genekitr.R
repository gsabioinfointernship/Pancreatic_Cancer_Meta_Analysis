library(genekitr)
library(tidyverse)
library(rio)

# Create folder for annotated results
dir.create("outputs/DESeq2/Annotated", showWarnings = FALSE)

# Load DESeq2 results
GSE130688 <- import("outputs/DESeq2/GSE130688.csv")

# Convert Gene_ID to character FIRST
GSE130688 <- GSE130688 %>%
  mutate(Gene_ID = as.character(Gene_ID))

# Get gene information including description and type
gene_info <- genInfo(GSE130688$Gene_ID, 
                     org = "hs", 
                     unique = TRUE, 
                     keepNA = FALSE)

# Merge all information - now all are character type
tab <- GSE130688 %>%
  left_join(gene_info, by = c("Gene_ID" = "input_id")) %>%
  dplyr::rename(Gene_Symbol = symbol,
                Gene_Description = gene_name,
                Gene_Biotype = gene_biotype) |>
  select(Gene_ID, Gene_Symbol, Gene_Description, Gene_Biotype, everything()) |> 
  filter(Gene_Biotype == "protein_coding")


# Export the annotated results
export(tab, "outputs/DESeq2/Annotated/GSE130688.csv")