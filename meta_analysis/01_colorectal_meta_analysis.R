source("meta_analysis/00_meta_analysis_functions.R")

cancer     <- "colorectal"
output_dir <- paste0("outputs/", cancer, "/meta_analysis/")

# Run meta analysis for colorectal cancer
run_meta_analysis(cancer, output_dir)
