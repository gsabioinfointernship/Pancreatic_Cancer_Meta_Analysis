library(survival)
library(survminer)
library(tidyverse)
library(rio)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Create output directories
dir.create("outputs/Survival", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/Survival", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/Survival/KM_curves", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Step 1: Load TCGA normalized expression data
# ============================================================================

normalized <- readRDS("data/TCGA/TCGA_PAAD_normalized.rds")
clinical <- import("data/TCGA/TCGA_PAAD_clinical.csv")

# ============================================================================
# Step 2: Filter for tumor samples with survival data
# ============================================================================

clinical_survival <- clinical |>
  filter(sample_type == "Tumor",!is.na(OS_days), OS_days > 0) |>
  select(sample_id, patient_id, OS_days, OS_status, age_years,
         gender, ajcc_pathologic_stage)

# ============================================================================
# Step 3: Load validated genes
# ============================================================================

validation <- import("outputs/TCGA/validation_results.csv")

# Step 4: Get top 50 validated genes
top_genes <- validation |>
  filter(validated == TRUE) |>
  arrange(padj) |>
  head(50) |>
  pull(Gene_ID)

# ============================================================================
# Step 5: Extract expression data for top genes
# ============================================================================

expression_data <- normalized[top_genes, clinical_survival$sample_id] |>
  t() |>
  as.data.frame() |>
  rownames_to_column("sample_id")

# ============================================================================
# Step 6: Combine expression with clinical data
# ============================================================================

survival_data <- clinical_survival |>
  left_join(expression_data, by = "sample_id")

# ============================================================================
# Step 7: Kaplan-Meier analysis for each gene
# ============================================================================

km_results <- tibble()

for (gene in top_genes) {
  if (gene %in% colnames(survival_data)) {

    # Calculate median expression
    median_expr <- median(survival_data[[gene]], na.rm = TRUE)

    # Create groups
    survival_data$group <- ifelse(survival_data[[gene]] >= median_expr, "High", "Low")

    # Fit survival model
    fit <- survfit(Surv(OS_days, OS_status) ~ group, data = survival_data)

    # Log-rank test
    diff <- survdiff(Surv(OS_days, OS_status) ~ group, data = survival_data)
    pval <- 1 - pchisq(diff$chisq, 1)

    # Cox model for HR
    cox_model <- coxph(Surv(OS_days, OS_status) ~ group, data = survival_data)
    hr <- exp(coef(cox_model))
    ci <- exp(confint(cox_model))

    # Store results
    km_results <- bind_rows(km_results, tibble(
      Gene_ID = gene,
      log_rank_p = pval,
      HR = hr,
      HR_lower = ci[1],
      HR_upper = ci[2],
      cox_p = summary(cox_model)$coefficients[5]
    ))
  }
}

# ============================================================================
# Step 8: Adjust p-values for multiple testing
# ============================================================================

km_results <- km_results |>
  mutate(
    FDR = p.adjust(log_rank_p, method = "BH"),
    significant = FDR < 0.05
  ) |>
  arrange(log_rank_p)

# Step 9: Export KM results
export(km_results, "outputs/Survival/KM_results.csv")

# ============================================================================
# Step 10: Generate KM plots for top 20 genes
# ============================================================================

top_20_survival <- km_results |>
  arrange(log_rank_p) |>
  head(20)

for (i in 1:nrow(top_20_survival)) {
  gene <- top_20_survival$Gene_ID[i]

  if (gene %in% colnames(survival_data)) {

    # Calculate median
    median_expr <- median(survival_data[[gene]], na.rm = TRUE)
    survival_data$group <- ifelse(survival_data[[gene]] >= median_expr, "High", "Low")

    # Fit survival
    fit <- survfit(Surv(OS_days, OS_status) ~ group, data = survival_data)

    # Create KM plot
    p <- ggsurvplot(
      fit = fit,
      data = survival_data,
      pval = TRUE,
      risk.table = TRUE,
      conf.int = FALSE,
      xlab = "Time (days)",
      ylab = "Overall Survival Probability",
      title = paste0(gene, " Expression"),
      legend.title = "Expression",
      legend.labs = c("High", "Low"),
      palette = c("#de2d26", "#2171b5"),
      ggtheme = theme_bw()
    )

    # Save plot
    ggsave(
      filename = paste0("figures/Survival/KM_curves/", gene, "_KM.png"),
      plot = print(p),
      width = 8,
      height = 8,
      dpi = 300
    )
  }
}

# ============================================================================
# Step 11: Univariate Cox regression for all genes
# ============================================================================

cox_univariate <- tibble()

for (gene in top_genes) {
  if (gene %in% colnames(survival_data)) {

    # Fit Cox model
    formula_cox <- as.formula(paste("Surv(OS_days, OS_status) ~", gene))
    cox_model <- coxph(formula_cox, data = survival_data)
    summary_cox <- summary(cox_model)

    # Extract results
    cox_univariate <- bind_rows(cox_univariate, tibble(
      Gene_ID = gene,
      HR = summary_cox$conf.int[1],
      HR_lower = summary_cox$conf.int[3],
      HR_upper = summary_cox$conf.int[4],
      p_value = summary_cox$coefficients[5],
      concordance = summary_cox$concordance[1]
    ))
  }
}

# Step 12: Adjust p-values
cox_univariate <- cox_univariate |>
  mutate(FDR = p.adjust(p_value, method = "BH")) |>
  arrange(p_value)

# Step 13: Export univariate Cox results
export(cox_univariate, "outputs/Survival/Cox_univariate.csv")

# ============================================================================
# Step 14: Multivariate Cox regression
# ============================================================================

# Get top prognostic genes
top_prognostic <- cox_univariate |>
  filter(FDR < 0.05) |>
  head(10) |>
  pull(Gene_ID)

if (length(top_prognostic) >= 3) {

  # Filter complete cases
  survival_data_complete <- survival_data |>
    filter(!is.na(age_years), !is.na(gender), !is.na(ajcc_pathologic_stage))

  if (nrow(survival_data_complete) > 50) {

    # Build formula
    formula_multi <- as.formula(paste(
      "Surv(OS_days, OS_status) ~",
      paste(top_prognostic, collapse = " + "),
      "+ age_years + gender + ajcc_pathologic_stage"
    ))

    # Fit multivariate Cox model
    cox_multivariate <- coxph(formula_multi, data = survival_data_complete)
    summary_multi <- summary(cox_multivariate)

    # Extract results
    multi_results <- as.data.frame(summary_multi$conf.int) |>
      rownames_to_column("Variable") |>
      mutate(p_value = summary_multi$coefficients[, 5]) |>
      select(Variable, `exp(coef)`, `lower .95`, `upper .95`, p_value) |>
      rename(HR = `exp(coef)`, HR_lower = `lower .95`, HR_upper = `upper .95`)

    # Export multivariate results
    export(multi_results, "outputs/Survival/Cox_multivariate.csv")
  }
}

# ============================================================================
# Step 15: Risk score development
# ============================================================================

if (length(top_prognostic) >= 3) {

  # Build risk model formula
  risk_formula <- as.formula(paste(
    "Surv(OS_days, OS_status) ~",
    paste(top_prognostic, collapse = " + ")
  ))

  # Fit risk model
  risk_model <- coxph(risk_formula, data = survival_data)

  # Calculate risk scores
  survival_data$risk_score <- predict(risk_model, type = "lp")

  # Create risk groups
  median_risk <- median(survival_data$risk_score)
  survival_data$risk_group <- ifelse(survival_data$risk_score >= median_risk,
                                     "High Risk", "Low Risk")

  # Fit survival by risk group
  fit_risk <- survfit(Surv(OS_days, OS_status) ~ risk_group, data = survival_data)

  # Create risk KM plot
  p_risk <- ggsurvplot(
    fit = fit_risk,
    data = survival_data,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = TRUE,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    title = "Risk Score Stratification",
    legend.title = "Risk Group",
    legend.labs = c("High Risk", "Low Risk"),
    palette = c("#de2d26", "#2171b5"),
    ggtheme = theme_bw()
  )

  # Save risk plot
  ggsave(
    filename = "figures/Survival/risk_score_KM.png",
    plot = print(p_risk),
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ============================================================================
# Step 16: Forest plot for top prognostic genes
# ============================================================================

forest_data <- cox_univariate |>
  filter(FDR < 0.05) |>
  head(30) |>
  mutate(
    Gene_ID = factor(Gene_ID, levels = rev(Gene_ID)),
    significant = p_value < 0.05
  )

if (nrow(forest_data) > 0) {

  p_forest <- ggplot(forest_data, aes(x = HR, y = Gene_ID)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.3) +
    geom_point(aes(color = significant), size = 3) +
    scale_color_manual(values = c("TRUE" = "#de2d26", "FALSE" = "gray60")) +
    scale_x_continuous(trans = "log10") +
    labs(
      title = "Forest Plot: Univariate Cox Regression",
      x = "Hazard Ratio (95% CI, log scale)",
      y = "Gene"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  # Save forest plot
  ggsave(
    filename = "figures/Survival/forest_plot.png",
    plot = p_forest,
    width = 10,
    height = 8,
    dpi = 300
  )
}
