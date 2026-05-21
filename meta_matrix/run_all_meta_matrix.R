# Run all meta_matrix scripts sequentially
scripts <- c(
  "meta_matrix/01_colorectal_meta_matrix.R",
  "meta_matrix/02_esophagus_meta_matrix.R",
  "meta_matrix/03_kidney_meta_matrix.R",
  "meta_matrix/04_liver_meta_matrix.R",
  "meta_matrix/05_pancreatic_meta_matrix.R"
)

for (script in scripts) {
  message("Running: ", script)
  source(script)
}

message("\nAll meta_matrix scripts complete.")
