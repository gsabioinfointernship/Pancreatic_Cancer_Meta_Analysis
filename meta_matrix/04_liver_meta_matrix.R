library(tidyverse)
library(rio)

cancer      <- "liver"
input_dir   <- paste0("outputs/", cancer, "/annotated/")
output_dir  <- paste0("outputs/", cancer, "/meta_matrix/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

geo_ids <- c(
  "GSE202853", "GSE207435", "GSE208036", "GSE216613",
  "GSE244826", "GSE263786", "GSE272035"
)

# Build named list and merged data frame
dataset_list <- list()
merged_list  <- list()

for (geo_id in geo_ids) {
  f <- paste0(input_dir, geo_id, ".csv")
  if (!file.exists(f)) {
    message("Skipping ", geo_id, ": annotated file not found")
    next
  }
  df <- import(f) |>
    mutate(GEO_ID = geo_id, Gene_ID = as.character(Gene_ID)) |>
    select(GEO_ID, everything())

  dataset_list[[geo_id]] <- df
  merged_list[[geo_id]]  <- df
}

# Format 1: merged long-format CSV
merged_df <- bind_rows(merged_list)
export(merged_df, paste0(output_dir, cancer, "_meta_matrix_merged.csv"))
message("Saved: ", output_dir, cancer, "_meta_matrix_merged.csv")

# Format 2: named list RDS
saveRDS(dataset_list, paste0(output_dir, cancer, "_meta_matrix_list.rds"))
message("Saved: ", output_dir, cancer, "_meta_matrix_list.rds")
