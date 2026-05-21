source("meta_analysis/00_meta_analysis_functions.R")

cancer     <- "liver"
output_dir <- paste0("outputs/", cancer, "/meta_analysis/")

message("Starting meta-analysis: ", cancer)
run_meta_analysis(cancer, output_dir)
message("Complete: ", cancer)
