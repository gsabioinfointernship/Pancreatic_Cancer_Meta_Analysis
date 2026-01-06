library(MetaVolcanoR)
library(tidyverse)
library(rio)

# Create output directories
dir.create("outputs/MetaVolcanoR", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/Meta-Analysis", showWarnings = FALSE, recursive = TRUE)

# Step 1: Define input directory
input_folder <- "outputs/DESeq2/Annotated/"

# Step 2: Get list of annotated DESeq2 result files
results_files <- list.files(
  path = input_folder,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Step 3: Extract GEO accession names
geo_accession <- tools::file_path_sans_ext(basename(results_files))

# Step 4: Read all files into a named list
studies <- import_list(results_files)
names(studies) <- geo_accession

# Step 5: Run random effects meta-analysis
meta_degs_rem <- rem_mv(
  diffexp = studies,
  pcriteria = "padj",
  foldchangecol = "log2FoldChange",
  genenamecol = "Gene_Symbol",
  geneidcol = "Gene_ID",
  collaps = TRUE,
  vcol = "lfcSE",
  cvar = FALSE,
  metathr = 0.05,
  jobname = "MetaVolcano",
  outputfolder = "figures/Meta-Analysis/",
  draw = "PDF",
  ncores = 6
)

# Step 6: Extract meta-analysis results
meta_results <- meta_degs_rem@metaresult

# Step 7: Export complete meta-analysis results
export(meta_results, "outputs/MetaVolcanoR/REM.csv")

# Step 8: Filter for significant genes
significant_degs <- meta_results |>
  filter(
    randomP < 0.05,
    abs(randomSummary) >= 1,
    abs(signcon) >= 2
  )

# Step 9: Export significant genes
export(significant_degs, "outputs/MetaVolcanoR/significant_genes.csv")
