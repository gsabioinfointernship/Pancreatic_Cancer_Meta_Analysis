library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)
library(rio)

# Create output directory
dir.create("data/TCGA", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Step 1: Query TCGA-PAAD RNA-seq data
# ============================================================================

query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# ============================================================================
# Step 2: Download TCGA data
# ============================================================================

GDCdownload(query)

# ============================================================================
# Step 3: Prepare data
# ============================================================================

data <- GDCprepare(query)

# ============================================================================
# Step 4: Extract count matrix
# ============================================================================

counts <- assay(data, "unstranded")

# ============================================================================
# Step 5: Extract sample metadata
# ============================================================================

coldata <- colData(data)

# ============================================================================
# Step 6: Download clinical data
# ============================================================================

clinical_query <- GDCquery_clinic(
  project = "TCGA-PAAD",
  type = "clinical"
)

# ============================================================================
# Step 7: Process clinical data
# ============================================================================

clinical_data <- clinical_query |>
  select(
    submitter_id,
    age_at_diagnosis,
    gender,
    ajcc_pathologic_stage,
    ajcc_pathologic_t,
    ajcc_pathologic_n,
    ajcc_pathologic_m,
    vital_status,
    days_to_death,
    days_to_last_follow_up
  )

# Step 8: Calculate survival endpoints
clinical_data <- clinical_data |>
  mutate(
    OS_status = ifelse(vital_status == "Dead", 1, 0),
    OS_days = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death),
    OS_years = OS_days / 365.25,
    age_years = age_at_diagnosis / 365.25
  )

# ============================================================================
# Step 9: Create sample information table
# ============================================================================

sample_info <- coldata |>
  as.data.frame() |>
  mutate(
    sample_id = colnames(counts),
    patient_id = substr(sample_id, 1, 12),
    sample_type = ifelse(grepl("-11", sample_id), "Normal", "Tumor")
  )

# Step 10: Join with clinical data
sample_info <- sample_info |>
  left_join(clinical_data, by = c("patient_id" = "submitter_id"))

# ============================================================================
# Step 11: Save raw count data
# ============================================================================

saveRDS(counts, "data/TCGA/TCGA_PAAD_counts.rds")

# ============================================================================
# Step 12: Save clinical information
# ============================================================================

export(sample_info, "data/TCGA/TCGA_PAAD_clinical.csv")

# ============================================================================
# Step 13: Normalize counts using VST
# ============================================================================

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ sample_type
)

# Filter low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# Extract normalized counts
normalized_counts <- assay(vsd)

# ============================================================================
# Step 14: Save normalized data
# ============================================================================

saveRDS(normalized_counts, "data/TCGA/TCGA_PAAD_normalized.rds")
