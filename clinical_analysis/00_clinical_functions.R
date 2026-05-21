library(tidyverse)
library(rio)
library(ggplot2)
library(ggpubr)

# ── Shared Configuration
CANCERS <- c("colorectal", "esophagus", "kidney", "liver", "pancreatic")

# ── Helper: Load TCGA Clinical Metadata
load_tcga_clinical <- function(cancer) {
  path <- file.path("data/tcga_curated", cancer, paste0(cancer, "_meta.csv"))
  clinical <- if (file.exists(path)) read.csv(path) else NULL
  
  # Check if stage column exists
  stage_col <- if (!is.null(clinical)) grep("pathologic_stage|tumor_stage|ajcc_pathologic_stage", colnames(clinical), value = TRUE)[1] else NA
  
  if (is.null(clinical) || is.na(stage_col)) {
    message("    Stage data missing in curated file. Downloading from GDC...")
    library(TCGAbiolinks)
    
    # TCGA project config (matching 01_per_cancer_ml.R if possible, or common ones)
    TCGA_PROJECTS <- list(
      colorectal = c("TCGA-COAD", "TCGA-READ"),
      esophagus  = "TCGA-ESCA",
      kidney     = c("TCGA-KIRC", "TCGA-KIRP", "TCGA-KICH"),
      liver      = "TCGA-LIHC",
      pancreatic = "TCGA-PAAD"
    )
    
    projects <- TCGA_PROJECTS[[cancer]]
    if (is.null(projects)) return(clinical)
    
    gdc_clin <- tryCatch({
      lapply(projects, function(p) {
        df <- GDCquery_clinic(p, type = "clinical")
        # Convert any list columns to character (e.g., 'sites_of_involvement' in ESCA)
        df %>% mutate(across(where(is.list), ~ sapply(., paste, collapse = "; ")))
      }) |> 
        bind_rows() |>
        rename(sample_id = submitter_id) # submitting id is patient barcode in GDCquery_clinic
    }, error = function(e) {
      message("    Error downloading GDC clinical data: ", e$message)
      NULL
    })
    
    if (is.null(gdc_clin)) return(clinical)
    
    # If we have curated clinical, merge to keep existing info but add stage
    if (!is.null(clinical)) {
      # In curated meta, sample_id is often the full barcode (15-16 chars)
      # In GDCquery_clinic, submitter_id is the patient barcode (12 chars)
      # Let's match on the first 12 characters of sample_id
      clinical$patient_id_short <- substr(clinical$sample_id, 1, 12)
      clinical <- left_join(clinical, gdc_clin, by = c("patient_id_short" = "sample_id"))
    } else {
      clinical <- gdc_clin
    }
  }
  
  clinical
}

# ── Helper: Load ML Risk Scores
load_risk_scores <- function(cancer) {
  path <- file.path("outputs", cancer, "machine_learning", paste0(cancer, "_patient_risk_scores.csv"))
  if (!file.exists(path)) return(NULL)
  read.csv(path)
}

# ── Core: Correlate Risk with Clinical Stage
run_clinical_correlation <- function(risk_df, clinical_df) {
  if (is.null(risk_df) || is.null(clinical_df)) return(NULL)
  
  # Match samples
  # risk_df has sample_id (likely full barcode) and patient_id (12 chars)
  # clinical_df after download has ajcc_pathologic_stage
  if ("ajcc_pathologic_stage" %in% colnames(clinical_df)) {
    # If we downloaded from GDC, match on patient_id
    # Ensure patient_id exists in clinical_df (it might be in bcr_patient_barcode or sample_id from my load function)
    if (!"patient_id" %in% colnames(clinical_df)) {
       clinical_df$patient_id <- substr(clinical_df$sample_id, 1, 12)
    }
    merged <- inner_join(risk_df, clinical_df, by = "patient_id")
  } else {
    merged <- inner_join(risk_df, clinical_df, by = "sample_id")
  }
  
  if (nrow(merged) == 0) {
     message("    WARNING: Merging risk and clinical data resulted in 0 rows.")
     return(NULL)
  }

  # Find Stage column
  stage_col <- grep("pathologic_stage|tumor_stage|ajcc_pathologic_stage", colnames(merged), value = TRUE)[1]
  if (is.na(stage_col)) {
    message("    WARNING: No stage column found in clinical data.")
    return(NULL)
  }
  
  # Clean Stage column (simplify to I, II, III, IV)
  merged <- merged |>
    mutate(Stage = case_when(
      grepl("Stage IV", .data[[stage_col]], ignore.case = TRUE) ~ "IV",
      grepl("Stage III", .data[[stage_col]], ignore.case = TRUE) ~ "III",
      grepl("Stage II", .data[[stage_col]], ignore.case = TRUE) ~ "II",
      grepl("Stage I", .data[[stage_col]], ignore.case = TRUE) ~ "I",
      TRUE ~ NA_character_
    )) |>
    filter(!is.na(Stage)) |>
    mutate(Stage = factor(Stage, levels = c("I", "II", "III", "IV")))
    
  if (nrow(merged) < 5) {
     message("    WARNING: Insufficient samples after stage filtering (n=", nrow(merged), ")")
     return(NULL)
  }

  merged
}

# ── Plot: Risk Score vs Stage Boxplot
save_stage_boxplot <- function(data, cancer, path) {
  if (is.null(data) || nrow(data) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  
  # Stats
  stats <- compare_means(risk_score ~ Stage, data = data, method = "kruskal.test")
  pval_label <- sprintf("Kruskal-Wallis p = %.3g", stats$p)
  
  p <- ggplot(data, aes(x = Stage, y = risk_score, fill = Stage)) +
    geom_violin(alpha = 0.3, color = NA) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.1, alpha = 0.2, size = 1) +
    scale_fill_brewer(palette = "Reds") +
    labs(
      title = paste(tools::toTitleCase(cancer), "— Risk vs. Clinical Stage"),
      subtitle = pval_label,
      x = "Pathologic Stage", y = "LASSO-Cox Risk Score"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none")
    
  ggsave(path, plot = p, width = 6, height = 5, dpi = 600, bg = "white")
}
