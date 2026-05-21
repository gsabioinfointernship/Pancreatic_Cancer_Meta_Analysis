# Run all meta-analysis scripts sequentially
scripts <- c(
  "meta_analysis/01_colorectal_meta_analysis.R",
  "meta_analysis/02_esophagus_meta_analysis.R",
  "meta_analysis/03_kidney_meta_analysis.R",
  "meta_analysis/04_liver_meta_analysis.R",
  "meta_analysis/05_pancreatic_meta_analysis.R"
)

for (script in scripts) {
  message("Running: ", script)
  source(script)
}

message("\nAll meta-analysis scripts complete.")
