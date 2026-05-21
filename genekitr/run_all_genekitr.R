# Run all genekitr annotation scripts sequentially
scripts <- c(
  "genekitr/01_colorectal_genekitr.R",
  "genekitr/02_esophagus_genekitr.R",
  "genekitr/03_kidney_genekitr.R",
  "genekitr/04_liver_genekitr.R",
  "genekitr/05_pancreatic_genekitr.R"
)

for (script in scripts) {
  message("Running: ", script)
  source(script)
}

message("\nAll genekitr annotation scripts complete.")
