library(genekitr)
library(tidyverse)
library(rio)

# Create output directory
dir.create("outputs/esophagus/annotated", showWarnings = FALSE, recursive = TRUE)

# Define all GEO dataset IDs
geo_ids <- c(
  "GSE104926", "GSE111011", "GSE119436", "GSE126304",
  "GSE130078", "GSE149609", "GSE164158", "GSE167488",
  "GSE194116", "GSE207526", "GSE213565", "GSE234304",
  "GSE269447", "GSE29968", "GSE32424", "GSE94869"
)

# Process each dataset
for (geo_id in geo_ids) {

  message("Processing ", geo_id, " ...")

  # Step 1: Load DESeq2 results
  input_file <- paste0("outputs/esophagus/DESeq2/", geo_id, ".csv")
  if (!file.exists(input_file)) {
    message("  Skipping ", geo_id, ": DESeq2 output not found")
    next
  }
  deseq_results <- import(input_file)

  # Step 2: Convert Gene_ID to character
  deseq_results <- deseq_results |>
    mutate(Gene_ID = as.character(Gene_ID))

  # Step 3: Get gene information from genekitr
  gene_info <- tryCatch(
    genInfo(id = deseq_results$Gene_ID, org = "hs", unique = TRUE, keepNA = FALSE),
    error = function(e) {
      message("  genInfo failed for ", geo_id, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(gene_info)) next

  # Step 4: Join gene information with DESeq2 results
  annotated <- deseq_results |>
    left_join(gene_info, by = c("Gene_ID" = "input_id"))

  # Step 5: Rename columns if present
  colnames(annotated)[colnames(annotated) == "symbol"]       <- "Gene_Symbol"
  colnames(annotated)[colnames(annotated) == "gene_name"]    <- "Gene_Description"
  colnames(annotated)[colnames(annotated) == "gene_biotype"] <- "Gene_Biotype"

  # Step 6: Keep only relevant columns
  annotated <- annotated |>
    select(any_of(c("Gene_ID", "Gene_Symbol", "Gene_Description", "Gene_Biotype",
                    "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")))

  # Step 7: Filter for protein-coding genes only
  if ("Gene_Biotype" %in% colnames(annotated)) {
    annotated <- annotated |>
      filter(Gene_Biotype == "protein_coding")
  }

  # Step 8: Export annotated results
  output_file <- paste0("outputs/esophagus/annotated/", geo_id, ".csv")
  export(annotated, output_file)
  message("  Done: ", output_file)
}
