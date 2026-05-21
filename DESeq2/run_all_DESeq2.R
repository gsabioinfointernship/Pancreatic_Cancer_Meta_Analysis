# Run all DESeq2 scripts sequentially
scripts <- c(
  "DESeq2/01_Colorectal_DESeq2.R",
  "DESeq2/02_Esophagus_DESeq2.R",
  "DESeq2/03_Kidney_DESeq2.R",
  "DESeq2/04_Liver_DESeq2.R",
  "DESeq2/05_Pancreatic_DESeq2.R"
)

for (script in scripts) {
  message("Running: ", script)
  source(script)
}

message("\nAll DESeq2 scripts complete.")
