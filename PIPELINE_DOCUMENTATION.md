# Pan-Cancer Meta-Analysis Pipeline — Documentation

> **Project:** Pan-Cancer Transcriptomic Meta-Analysis (5 Cancers)
> **Cancers:** Colorectal, Esophagus, Kidney, Liver, Pancreatic
> **Platform:** R (Bioconductor)

---

## Table of Contents

1. [Project Structure](#1-project-structure)
2. [Data Sources](#2-data-sources)
3. [Step 1 — DESeq2 Differential Expression](#3-step-1--deseq2-differential-expression)
4. [Step 2 — Meta Matrix Preparation](#4-step-2--meta-matrix-preparation)
5. [Step 3 — Meta-Analysis (Random Effects Model)](#5-step-3--meta-analysis-random-effects-model)
6. [Step 4 — Functional Enrichment](#6-step-4--functional-enrichment)
7. [Step 5 — Network Analysis](#7-step-5--network-analysis)
8. [Step 6 — Drug Target Analysis](#8-step-6--drug-target-analysis)
9. [Step 7 — Survival Analysis](#9-step-7--survival-analysis)
10. [Step 8 — Validation (TCGA)](#10-step-8--validation-tcga)
11. [Step 9 — Machine Learning](#11-step-9--machine-learning)
12. [Step 10 — Immune Microenvironment Analysis](#12-step-10--immune-microenvironment-analysis)
13. [Step 11 — Master Regulator Analysis](#13-step-11--master-regulator-analysis)
14. [Step 12 — Clinicopathological Correlation](#14-step-12--clinicopathological-correlation)
15. [Key Parameters & Thresholds](#15-key-parameters--thresholds)
16. [Output File Reference](#16-output-file-reference)
17. [How to Re-Run](#17-how-to-re-run)

---

## 1. Project Structure

```
Pan_Cancer_Meta/
├── data/
│   ├── {cancer}/raw_counts/       ← GEO raw count matrices (TSV/CSV)
│   ├── {cancer}/metadata/         ← sample metadata per dataset
│   └── tcga_curated/              ← cached TCGA downloads (curatedTCGAData)
│
├── outputs/
│   └── {cancer}/
│       ├── DESeq2/                ← per-dataset DESeq2 results
│       ├── annotated/             ← DESeq2 results with gene symbols
│       ├── meta_matrix/           ← combined input for meta-analysis
│       ├── meta_analysis/         ← pooled effect sizes, p-values
│       ├── enrichment/            ← ORA/GSEA tables
│       ├── network/               ← PPI metrics, hub genes, modules
│       ├── drug_targets/          ← DGIdb interactions, priority table
│       ├── survival/              ← Cox results, prognostic genes
│       ├── validation/            ← TCGA limma validation results
│       ├── machine_learning/      ← RF/LASSO models and metrics
│       ├── immune_analysis/       ← ssGSEA scores and correlation tables
│       ├── regulatory_analysis/   ← TF enrichment and GSEA
│       ├── clinical_analysis/     ← Risk vs. Stage/Grade correlations
│       └── figures/               ← all plots (subfolders per analysis)
│           ├── enrichment/
│           ├── network/
│           ├── drug_targets/
│           ├── survival/
│           ├── validation/
│           ├── machine_learning/  ← ROC, KM plots, importance
│           ├── immune_analysis/   ← Infiltration landscapes, dotplots
│           ├── regulatory_analysis/ ← TF dotplots, heatmaps
│           └── clinical_analysis/   ← Stage boxplots, grid plots
│
├── DESeq2/                        ← per-cancer DESeq2 scripts
├── meta_matrix/                   ← per-cancer meta-matrix scripts
├── meta_analysis/                 ← meta-analysis scripts
├── functional_enrichment/         ← enrichment scripts
├── network_analysis/              ← network scripts
├── drug_target_analysis/          ← drug target scripts
├── survival_analysis/             ← survival analysis scripts
├── validation/                    ← TCGA validation scripts
├── machine_learning/              ← ML diagnostic/prognostic scripts
├── immune_analysis/               ← Immune infiltration scripts (ssGSEA)
├── regulatory_analysis/           ← Transcription Factor analysis
├── clinical_analysis/             ← Clinicopathological analysis
└── scripts/                       ← legacy/utility scripts
```

---

## 13. Step 11 — Master Regulator Analysis

**Scripts:** `regulatory_analysis/00_regulatory_functions.R`, `01_per_cancer_regulatory.R`, `02_pancancer_regulatory.R`
**Runner:** `scripts/run_advanced_deepdives.R` (Step 1)

### What it does
- Identifies the **Transcription Factors (TFs)** that act as "Master Regulators" of the meta-analysis signatures.
- Uses the **MSigDB C3:TFT (GTRD)** database of TF-target interactions.
- Performs **ORA** on up/down DEGs and **GSEA** on the full ranked list.
- Outputs: TF enrichment tables, dotplots of top regulators, and pan-cancer consistent TF heatmaps.

### Impact
Moves the study from "descriptive" to "regulatory," identifying the upstream switches driving the cancer transcriptome.

---

## 14. Step 12 — Clinicopathological Correlation

**Scripts:** `clinical_analysis/00_clinical_functions.R`, `01_per_cancer_clinical.R`, `02_pancancer_clinical.R`
**Runner:** `scripts/run_advanced_deepdives.R` (Step 2)

### What it does
- Directly correlates the **ML Risk Scores** (LASSO-Cox) with patient clinical features in TCGA.
- Focuses on **Tumor Stage (I-IV)** and progression.
- Performs **Kruskal-Wallis** statistical tests to prove risk score increases with disease severity.
- Outputs: Merged risk-clinical tables, Stage boxplots, and pan-cancer grid plots.

### Impact
Provides clinical proof that the identified molecular signature tracks with actual disease progression in patients.

---

## 15. Key Parameters & Thresholds

| Parameter | Value | Used In |
|-----------|-------|---------|
| `PADJ_CUTOFF` | 0.05 | All steps |
| `LFC_CUTOFF` | 1 (log2FC) | Meta-analysis, enrichment, network, drug |
| `STRING_SCORE` | 700 / 1000 | Network analysis |
| `HUB_QUANTILE` | 0.90 (top 10%) | Network analysis |
| `MIN_MODULE_SIZE` | 5 nodes | Network analysis |
| `MIN_CANCERS_PC` | 3 | Pan-cancer meta, network, drug |
| `HUB_N` | 100 | Survival & ML — max candidate genes |
| `COX_PADJ` | 0.05 | Survival — prognostic significance |
| `MIN_EVENTS` | 10-15 | Survival & ML — min deaths required |
| `RF_TREES` | 500 | ML — random forest depth |
| `TRAIN_FRAC` | 0.8 | ML — train/test split |
| `GLMNET_MAXIT` | 1,000,000 | ML — convergence limit |

---

## 16. Output File Reference
*(See previous documentation)*

---

## 17. How to Re-Run

### Full pipeline (all steps)
```r
# Step 1-3: Expression & Meta-analysis
source("DESeq2/run_all_DESeq2.R")
source("meta_matrix/run_all_meta_matrix.R")
source("meta_analysis/run_all_meta_analysis.R")

# Step 4-6: Downstream Biology
source("functional_enrichment/run_all_enrichment.R")
source("network_analysis/run_all_network.R")
source("drug_target_analysis/run_all_drug.R")

# Step 7-10: Clinical Validation, ML & Immune
source("survival_analysis/run_all_survival.R")
source("validation/run_all_validation.R")
source("machine_learning/run_all_ml.R")
source("immune_analysis/run_all_immune.R")

# Step 11-12: Advanced Deep-Dives
source("scripts/run_advanced_deepdives.R")
```

---

*Documentation generated for Pan_Cancer_Meta project. Last updated: 2026-05-05.*
