source("meta_analysis/00_meta_analysis_functions.R")

cancer     <- "esophagus"
output_dir <- paste0("outputs/", cancer, "/meta_analysis/")

message("Starting meta-analysis: ", cancer)
run_meta_analysis(cancer, output_dir)
message("Complete: ", cancer)
