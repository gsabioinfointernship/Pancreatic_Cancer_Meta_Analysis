setwd("/Users/jubayer/Projects/GSA/PAAD_Meta")

cat("\n================================================================================\n")
cat("PAAD Meta-Analysis: Complete Workflow Execution\n")
cat("================================================================================\n\n")

start_time_total <- Sys.time()

run_script <- function(script_name, description) {
  cat("\n--------------------------------------------------------------------------------\n")
  cat("RUNNING:", description, "\n")
  cat("Script:", script_name, "\n")
  cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("--------------------------------------------------------------------------------\n")

  start_time <- Sys.time()

  tryCatch({
    source(script_name)
    end_time <- Sys.time()
    elapsed <- round(difftime(end_time, start_time, units = "mins"), 2)

    cat("\n✓ COMPLETED:", description, "\n")
    cat("  Runtime:", elapsed, "minutes\n")

    return(list(success = TRUE, time = elapsed))

  }, error = function(e) {
    cat("\n✗ FAILED:", description, "\n")
    cat("  Error:", conditionMessage(e), "\n")

    return(list(success = FALSE, time = NA, error = conditionMessage(e)))
  })
}

results <- list()

cat("\n🔹 PHASE 1: DATA PREPARATION\n")
results[["01_genekitr"]] <- run_script(
  "scripts/01_run_genekitr_all.R",
  "Batch Gene Annotation (CRITICAL PREREQUISITE)"
)

if(!results[["01_genekitr"]]$success) {
  stop("Script 01 failed. Cannot proceed without annotated data.")
}

cat("\n🔹 PHASE 2: FUNCTIONAL ENRICHMENT\n")
results[["02_enrichment"]] <- run_script(
  "scripts/02_functional_enrichment.R",
  "Functional Enrichment Analysis"
)

cat("\n🔹 PHASE 3: TCGA DATA & VALIDATION\n")
results[["03_tcga_download"]] <- run_script(
  "scripts/03_TCGA_data_download.R",
  "TCGA-PAAD Data Download"
)

if(results[["03_tcga_download"]]$success) {
  results[["04_tcga_validation"]] <- run_script(
    "scripts/04_TCGA_validation.R",
    "TCGA Validation Analysis"
  )
} else {
  cat("\n⚠ Skipping TCGA validation (download failed)\n")
  results[["04_tcga_validation"]] <- list(success = FALSE, time = NA, error = "Skipped")
}

cat("\n🔹 PHASE 4: SURVIVAL ANALYSIS\n")
if(results[["03_tcga_download"]]$success) {
  results[["05_survival"]] <- run_script(
    "scripts/05_survival_analysis.R",
    "Survival & Prognostic Analysis"
  )
} else {
  cat("\n⚠ Skipping survival analysis (TCGA data not available)\n")
  results[["05_survival"]] <- list(success = FALSE, time = NA, error = "Skipped")
}

cat("\n🔹 PHASE 5: MACHINE LEARNING\n")
if(results[["03_tcga_download"]]$success) {
  results[["06_ml"]] <- run_script(
    "scripts/06_machine_learning.R",
    "Machine Learning Classification"
  )
} else {
  cat("\n⚠ Skipping machine learning (TCGA data not available)\n")
  results[["06_ml"]] <- list(success = FALSE, time = NA, error = "Skipped")
}

cat("\n🔹 PHASE 6: NETWORK & DRUG ANALYSIS\n")
results[["07_network"]] <- run_script(
  "scripts/07_network_analysis.R",
  "PPI Network Analysis"
)

results[["08_drug"]] <- run_script(
  "scripts/08_drug_target_analysis.R",
  "Drug Target Identification"
)

cat("\n🔹 PHASE 7: CLINICAL NOMOGRAM (OPTIONAL)\n")
if(results[["05_survival"]]$success) {
  results[["09_nomogram"]] <- run_script(
    "scripts/09_clinical_nomogram.R",
    "Clinical Prognostic Nomogram"
  )
} else {
  cat("\n⚠ Skipping nomogram (survival analysis not completed)\n")
  results[["09_nomogram"]] <- list(success = FALSE, time = NA, error = "Skipped")
}

cat("\n🔹 PHASE 8: PUBLICATION OUTPUTS\n")
results[["10_pub_figs"]] <- run_script(
  "scripts/10_publication_figures.R",
  "Main Publication Figures"
)

results[["11_supp_figs"]] <- run_script(
  "scripts/11_supplementary_figures.R",
  "Supplementary Figures"
)

results[["12_supp_tables"]] <- run_script(
  "scripts/12_supplementary_tables.R",
  "Supplementary Tables"
)

end_time_total <- Sys.time()
total_elapsed <- round(difftime(end_time_total, start_time_total, units = "mins"), 2)

cat("\n================================================================================\n")
cat("ANALYSIS WORKFLOW SUMMARY\n")
cat("================================================================================\n\n")

summary_df <- data.frame(
  Phase = c(
    "01. Gene Annotation",
    "02. Enrichment Analysis",
    "03. TCGA Download",
    "04. TCGA Validation",
    "05. Survival Analysis",
    "06. Machine Learning",
    "07. Network Analysis",
    "08. Drug Targets",
    "09. Clinical Nomogram",
    "10. Publication Figures",
    "11. Supplementary Figures",
    "12. Supplementary Tables"
  ),
  Status = sapply(results, function(x) ifelse(x$success, "✓ Success", "✗ Failed")),
  Runtime_Minutes = sapply(results, function(x) ifelse(is.na(x$time), "-", as.character(x$time)))
)

print(summary_df, row.names = FALSE)

cat("\n--------------------------------------------------------------------------------\n")
cat("Total Runtime:", total_elapsed, "minutes (", round(total_elapsed/60, 2), "hours)\n")
cat("Completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

success_count <- sum(sapply(results, function(x) x$success))
total_count <- length(results)

cat("\nSuccessful scripts:", success_count, "/", total_count, "\n")

if(success_count == total_count) {
  cat("\n🎉 ALL ANALYSES COMPLETED SUCCESSFULLY!\n")
  cat("\nNext steps:\n")
  cat("  1. Review all outputs in outputs/ and figures/ directories\n")
  cat("  2. Check quality control metrics in each output folder\n")
  cat("  3. Proceed with manuscript preparation\n")
  cat("  4. See RESEARCH_PLAN.md for detailed interpretation\n")
} else {
  cat("\n⚠ SOME ANALYSES FAILED OR WERE SKIPPED\n")
  cat("\nFailed/Skipped scripts:\n")
  for(name in names(results)) {
    if(!results[[name]]$success) {
      cat("  -", name, "\n")
      if(!is.null(results[[name]]$error)) {
        cat("    Error:", results[[name]]$error, "\n")
      }
    }
  }
  cat("\nRefer to ANALYSIS_WORKFLOW.md for troubleshooting\n")
}

cat("\n================================================================================\n")
