# Meta-Analysis of RNA-seq data using MetaVolcanoR
# Author: Muntasim Fuad

# Load required pacakges
library(MetaVolcanoR) 
library(tidyverse)
library(rio)

# Define input directory
input_folder <- "outputs/DESeq2/Annotated/"

# Generate list of CSV file paths and extract project names
results <- list.files(path = input_folder, pattern = "\\.csv$", full.names = TRUE)
geo_accession <- tools::file_path_sans_ext(basename(results))

# Read all files into a named list using rio
studies <- import_list(results)
names(studies) <- geo_accession

# ------------------------------------------------------------------------------
# Meta-Analysis

# Create a output folder for figures
dir.create("figures/Meta-Analysis/", showWarnings = FALSE)

# Random Effect Model
meta_degs_rem <- rem_mv(diffexp= studies,
                        pcriteria='padj',
                        foldchangecol= "log2FoldChange",
                        genenamecol= "Gene_Symbol",
                        geneidcol= "Gene_ID",
                        collaps= TRUE,
                        vcol= "lfcse", 
                        cvar=FALSE,
                        metathr= 0.05,
                        jobname= "MetaVolcano",
                        outputfolder= "figures/Meta-Analysis/", 
                        draw= 'PDF',
                        ncores= 6)

meta_results <- meta_degs_rem@metaresult

# Export the results
export(meta_results, "outputs/MetaVolcanoR/REM.csv")

# Filter results 
# Filter results 
sign_degs <- meta_results |> 
  filter(randomP < 0.05, 
         abs(randomSummary) >= 1,
         abs(signcon) >= 2
  )
        
# Export the results 
export(sign_degs, "outputs/MetaVolcanoR/significant_genes.csv")

