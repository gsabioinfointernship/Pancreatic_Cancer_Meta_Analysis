# PAAD Meta-Analysis: Complete Analysis Workflow

## Overview
This workflow executes a comprehensive RNA-seq meta-analysis for pancreatic adenocarcinoma (PAAD) biomarker discovery, validation, and therapeutic target identification.

## Prerequisites

### R Packages Installation
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "sva", "clusterProfiler", "enrichplot",
  "org.Hs.eg.db", "DOSE", "ReactomePA", "pathview",
  "TCGAbiolinks", "SummarizedExperiment", "ComplexHeatmap",
  "STRINGdb", "genekitr"
))

install.packages(c(
  "tidyverse", "rio", "MetaVolcanoR", "survival", "survminer",
  "caret", "randomForest", "e1071", "glmnet", "xgboost",
  "pROC", "ROCR", "igraph", "ggraph", "visNetwork", "rms",
  "patchwork", "cowplot", "circlize", "pheatmap", "ggpubr",
  "msigdbr", "rDGIdb", "writexl"
))
```

---

## Execution Order

### Phase 0: Individual Differential Expression Analysis (REQUIRED FIRST)

**Script 00a: DESeq2 Analysis for All Datasets**
```r
source("scripts/00_DESeq2_all.R")
```
- Processes all 7 GEO datasets individually
- DESeq2 differential expression with SVA batch correction
- Generates PCA plots (pre/post correction) for quality control
- **Creates**: `outputs/DESeq2/` (7 CSV files), `figures/PCA/` (14 PNG files)
- **Runtime**: ~15-20 minutes
- **Critical**: Must run before meta-analysis

**Script 00b: Meta-Analysis with MetaVolcanoR**
```r
source("scripts/00_MetaVolcanoR.R")
```
- Combines results from all 7 datasets using random effects model
- Identifies consistently differentially expressed genes
- **Creates**: `outputs/MetaVolcanoR/` (REM.csv, significant_genes.csv)
- **Requires**: Script 00a and Script 01 must be completed first
- **Runtime**: ~5-10 minutes

---

### Phase 1: Data Preparation

**Script 01: Batch Gene Annotation**
```r
source("scripts/01_run_genekitr_all.R")
```
- Annotates all 7 DESeq2 results with gene symbols, descriptions, and biotypes
- Filters for protein-coding genes only
- **Creates**: `outputs/DESeq2/Annotated/` (7 CSV files)
- **Runtime**: ~5-10 minutes
- **Critical**: Required for meta-analysis (Script 00b) and all downstream analyses

---

### Phase 2: Functional Enrichment Analysis

**Script 02: Pathway Enrichment**
```r
source("scripts/02_functional_enrichment.R")
```
- GO, KEGG, Reactome pathway enrichment
- Gene Set Enrichment Analysis (GSEA)
- MSigDB Hallmark gene sets
- **Creates**: `outputs/Enrichment/`, `figures/Enrichment/`
- **Runtime**: ~10-15 minutes
- **Depends on**: Script 00b

---

### Phase 3: TCGA Data Acquisition & Validation

**Script 03: Download TCGA-PAAD Data**
```r
source("scripts/03_TCGA_data_download.R")
```
- Downloads TCGA-PAAD RNA-seq and clinical data
- Normalizes expression data
- **Creates**: `data/TCGA/` (counts, clinical, normalized RDS files)
- **Runtime**: ~20-30 minutes (depends on internet speed)
- **Note**: Requires TCGAbiolinks package

**Script 04: TCGA Validation Analysis**
```r
source("scripts/04_TCGA_validation.R")
```
- Validates 613 DEGs in TCGA-PAAD cohort
- Calculates correlation and concordance metrics
- Generates validation figures
- **Creates**: `outputs/TCGA/`, `figures/TCGA/`
- **Runtime**: ~5-10 minutes
- **Depends on**: Script 00b, Script 03

---

### Phase 4: Survival Analysis

**Script 05: Survival & Prognostic Analysis**
```r
source("scripts/05_survival_analysis.R")
```
- Kaplan-Meier survival curves for top 50 genes
- Univariate and multivariate Cox regression
- Risk score development
- **Creates**: `outputs/Survival/`, `figures/Survival/`
- **Runtime**: ~15-20 minutes
- **Depends on**: Script 03, Script 04

---

### Phase 5: Machine Learning Classification

**Script 06: ML Model Development**
```r
source("scripts/06_machine_learning.R")
```
- Trains Random Forest, SVM, Elastic Net, XGBoost models
- 10-fold cross-validation on GEO meta-cohort
- Independent validation on TCGA-PAAD
- **Creates**: `outputs/MachineLearning/`, `figures/MachineLearning/`
- **Runtime**: ~20-30 minutes
- **Depends on**: Script 00b, Script 03

---

### Phase 6: Network & Drug Analysis

**Script 07: PPI Network Analysis**
```r
source("scripts/07_network_analysis.R")
```
- Protein-protein interaction network construction
- Hub gene identification
- Network module detection
- **Creates**: `outputs/Network/`, `figures/Network/`
- **Runtime**: ~10-15 minutes
- **Depends on**: Script 00b

**Script 08: Drug Target Identification**
```r
source("scripts/08_drug_target_analysis.R")
```
- Drug-gene interaction database queries
- Druggability assessment
- Drug repurposing candidates
- **Creates**: `outputs/DrugTarget/`, `figures/DrugTarget/`
- **Runtime**: ~5-10 minutes
- **Depends on**: Script 00b, Script 07

---

### Phase 7: Clinical Nomogram (Optional)

**Script 09: Prognostic Nomogram**
```r
source("scripts/09_clinical_nomogram.R")
```
- Develops clinical prediction nomogram
- Calibration curves and risk stratification
- **Creates**: `outputs/Clinical/`, `figures/Clinical/`
- **Runtime**: ~5-10 minutes
- **Depends on**: Script 05
- **Note**: Requires sufficient prognostic genes (≥3 with FDR < 0.05)

---

### Phase 8: Publication Outputs

**Script 10: Main Publication Figures**
```r
source("scripts/10_publication_figures.R")
```
- Generates 8 main publication figures
- High-resolution PDFs suitable for journals
- **Creates**: `figures/Publication/`
- **Runtime**: ~5-10 minutes
- **Depends on**: All previous scripts

**Script 11: Supplementary Figures**
```r
source("scripts/11_supplementary_figures.R")
```
- Generates supplementary figures
- Individual study volcanos, extended enrichment plots, etc.
- **Creates**: `figures/Supplementary/`
- **Runtime**: ~5-10 minutes
- **Depends on**: All previous scripts

**Script 12: Supplementary Tables**
```r
source("scripts/12_supplementary_tables.R")
```
- Compiles all supplementary tables
- Creates Excel workbook with all tables
- **Creates**: `outputs/SupplementaryTables/`
- **Runtime**: ~2-5 minutes
- **Depends on**: All previous scripts

---

## Complete Workflow Execution

### Option 1: Run All Steps Sequentially

```r
setwd("/Users/jubayer/Projects/GSA/PAAD_Meta")

# Phase 0: Individual analysis and meta-analysis
source("scripts/00_DESeq2_all.R")        # Step 0a: REQUIRED FIRST
source("scripts/01_run_genekitr_all.R")  # Step 1: Annotate genes
source("scripts/00_MetaVolcanoR.R")      # Step 0b: Meta-analysis

# Phase 2: Enrichment
source("scripts/02_functional_enrichment.R")

# Phase 3: TCGA validation
source("scripts/03_TCGA_data_download.R")
source("scripts/04_TCGA_validation.R")

# Phase 4: Survival
source("scripts/05_survival_analysis.R")

# Phase 5: Machine learning
source("scripts/06_machine_learning.R")

# Phase 6: Network and drugs
source("scripts/07_network_analysis.R")
source("scripts/08_drug_target_analysis.R")

# Phase 7: Clinical nomogram
source("scripts/09_clinical_nomogram.R")

# Phase 8: Publication outputs
source("scripts/10_publication_figures.R")
source("scripts/11_supplementary_figures.R")
source("scripts/12_supplementary_tables.R")
```

**Total Runtime**: Approximately 2.5-3.5 hours

---

## Output Structure

```
PAAD_Meta/
├── outputs/
│   ├── DESeq2/                 # Individual dataset results
│   │   ├── GSE*.csv            # 7 files
│   │   └── Annotated/          # 7 annotated files
│   ├── MetaVolcanoR/           # Meta-analysis results
│   │   ├── REM.csv             # All genes
│   │   └── significant_genes.csv  # 613 DEGs
│   ├── Enrichment/             # GO, KEGG, Reactome results
│   ├── TCGA/                   # Validation results
│   ├── Survival/               # KM and Cox regression results
│   ├── MachineLearning/        # Model performance metrics
│   ├── Network/                # PPI network and hub genes
│   ├── DrugTarget/             # Druggable genes and interactions
│   ├── Clinical/               # Nomogram model
│   └── SupplementaryTables/    # All supplementary tables + Excel
│
├── figures/
│   ├── PCA/                    # QC plots for 7 datasets
│   ├── Meta-Analysis/          # Meta-volcano plot
│   ├── Enrichment/
│   ├── TCGA/
│   ├── Survival/
│   ├── MachineLearning/
│   ├── Network/
│   ├── DrugTarget/
│   ├── Clinical/
│   ├── Publication/            # 8 main figures
│   └── Supplementary/          # 11+ supplementary figures
│
└── data/
    └── TCGA/                   # TCGA-PAAD downloaded data
```

---

## Critical Workflow Steps

### Step 0a (DESeq2) is ABSOLUTELY REQUIRED
- Processes raw count data for all 7 datasets
- Performs batch correction using SVA
- Generates quality control PCA plots
- Without this, meta-analysis cannot proceed

### Step 1 (Annotation) is REQUIRED before Meta-Analysis
- Adds gene symbols and descriptions to DESeq2 results
- Required input for MetaVolcanoR (Script 00b)

### Correct Execution Order
1. **00_DESeq2_all.R** (individual analysis)
2. **01_run_genekitr_all.R** (gene annotation)
3. **00_MetaVolcanoR.R** (meta-analysis)
4. All other scripts (02-12) in order

---

## Expected Key Results

### Individual DESeq2 Analysis (Script 00a)
- 7 differential expression result files
- 14 PCA plots showing batch correction effectiveness

### Meta-Analysis (Script 00b)
- **613 significant DEGs** (padj < 0.05, |log2FC| ≥ 1, signcon ≥ 2)
- Top downregulated: CLCA1, APOA1, APOA4, TMPRSS15

### TCGA Validation
- Expected correlation: r > 0.5
- Direction concordance: >80%

### Survival Analysis
- Expected: 10-30 prognostic genes (FDR < 0.05)
- Multi-gene risk score for patient stratification

### Machine Learning
- Training AUC: 0.95-0.99
- Validation AUC: 0.90-0.95 (target: >0.90)

### Network Analysis
- 20-50 hub genes
- Multiple network modules

### Drug Targets
- 50-100 druggable genes
- 10-20 FDA-approved drugs (if rDGIdb available)

---

## Quality Control Checkpoints

After each phase, verify:

1. **After Script 00a**: Check `outputs/DESeq2/` has 7 files and `figures/PCA/` has 14 plots
2. **After Script 01**: Check `outputs/DESeq2/Annotated/` has 7 files
3. **After Script 00b**: Verify 613 significant genes in `outputs/MetaVolcanoR/significant_genes.csv`
4. **After Script 03**: Confirm TCGA download completed (~180-190 tumor samples)
5. **After Script 04**: Validation correlation r > 0.5, concordance > 70%
6. **After Script 05**: At least 5-10 prognostic genes identified
7. **After Script 06**: ML validation AUC > 0.85 (ideally >0.90)
8. **After Script 07**: Hub genes identified (>20 genes)
9. **After Script 08**: Druggable genes found (>30 genes)

---

## Troubleshooting

### Script 00a Issues
- **Low gene counts**: Check if raw count files are in correct format
- **SVA errors**: May indicate insufficient sample size or no batch effects

### Script 00b Issues
- **Meta-analysis fails**: Ensure Script 01 completed successfully (annotated files exist)
- **Few significant genes**: Check individual dataset quality

### Script 03 Takes Too Long
- **Issue**: TCGA download is slow
- **Solution**: Normal for large datasets, wait patiently or use faster internet

### Other Issues
- See original ANALYSIS_WORKFLOW.md for additional troubleshooting

---

## Citation

If you use this workflow, please cite:
- DESeq2: Love et al., Genome Biology 2014
- SVA: Leek et al., PLoS Genetics 2012
- MetaVolcanoR: Prada et al., GigaScience 2018
- clusterProfiler: Wu et al., Innovation 2021
- TCGAbiolinks: Colaprico et al., Nucleic Acids Research 2016

---

Last updated: January 2026
