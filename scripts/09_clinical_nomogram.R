library(rms)
library(survival)
library(tidyverse)
library(rio)

dir.create("outputs/Clinical", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/Clinical", showWarnings = FALSE, recursive = TRUE)

normalized <- readRDS("data/TCGA/TCGA_PAAD_normalized.rds")
clinical <- import("data/TCGA/TCGA_PAAD_clinical.csv")
cox_results <- import("outputs/Survival/Cox_univariate.csv")

top_prognostic <- cox_results %>%
  filter(FDR < 0.05) %>%
  arrange(p_value) %>%
  head(10) %>%
  pull(Gene_ID)

if(length(top_prognostic) >= 3) {

  clinical_survival <- clinical %>%
    filter(sample_type == "Tumor", !is.na(OS_days), OS_days > 0) %>%
    filter(!is.na(age_years), !is.na(ajcc_pathologic_stage)) %>%
    dplyr::select(sample_id, patient_id, OS_days, OS_status, age_years,
                  gender, ajcc_pathologic_stage)

  expression_data <- normalized[top_prognostic, clinical_survival$sample_id] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample_id")

  nomogram_data <- clinical_survival %>%
    left_join(expression_data, by = "sample_id") %>%
    mutate(
      OS_years = OS_days / 365.25,
      stage_simple = case_when(
        grepl("Stage I", ajcc_pathologic_stage) ~ "I",
        grepl("Stage II", ajcc_pathologic_stage) ~ "II",
        grepl("Stage III|Stage IV", ajcc_pathologic_stage) ~ "III-IV",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(stage_simple)) %>%
    na.omit()

  if(nrow(nomogram_data) >= 50) {

    dd <- datadist(nomogram_data)
    options(datadist = "dd")

    formula_genes <- paste(top_prognostic[1:min(5, length(top_prognostic))], collapse = " + ")
    formula_str <- paste("Surv(OS_years, OS_status) ~", formula_genes, "+ age_years + stage_simple")

    cox_model <- cph(as.formula(formula_str), data = nomogram_data, x = TRUE, y = TRUE, surv = TRUE)

    png("figures/Clinical/nomogram.png", width = 12, height = 8, units = "in", res = 300)
    nom <- nomogram(cox_model,
                    fun = list(function(x) survest(cox_model, times = 365.25, newdata = data.frame(x))$surv,
                               function(x) survest(cox_model, times = 365.25 * 3, newdata = data.frame(x))$surv,
                               function(x) survest(cox_model, times = 365.25 * 5, newdata = data.frame(x))$surv),
                    funlabel = c("1-Year Survival", "3-Year Survival", "5-Year Survival"))
    plot(nom)
    dev.off()

    cal_1yr <- calibrate(cox_model, u = 365.25, B = 100)
    cal_3yr <- calibrate(cox_model, u = 365.25 * 3, B = 100)
    cal_5yr <- calibrate(cox_model, u = 365.25 * 5, B = 100)

    png("figures/Clinical/calibration_curves.png", width = 12, height = 4, units = "in", res = 300)
    par(mfrow = c(1, 3))
    plot(cal_1yr, main = "1-Year Survival Calibration")
    plot(cal_3yr, main = "3-Year Survival Calibration")
    plot(cal_5yr, main = "5-Year Survival Calibration")
    dev.off()

    nomogram_data$predicted_risk <- predict(cox_model, type = "lp")

    risk_tertiles <- quantile(nomogram_data$predicted_risk, probs = c(1/3, 2/3))
    nomogram_data$risk_group <- cut(nomogram_data$predicted_risk,
                                     breaks = c(-Inf, risk_tertiles[1], risk_tertiles[2], Inf),
                                     labels = c("Low", "Medium", "High"))

    fit_risk <- survfit(Surv(OS_years, OS_status) ~ risk_group, data = nomogram_data)

    png("figures/Clinical/risk_stratification.png", width = 10, height = 8, units = "in", res = 300)
    plot(fit_risk,
         col = c("#2171b5", "#fdae61", "#de2d26"),
         lwd = 2,
         xlab = "Time (years)",
         ylab = "Overall Survival Probability",
         main = "Risk Stratification by Nomogram")
    legend("topright",
           legend = c("Low Risk", "Medium Risk", "High Risk"),
           col = c("#2171b5", "#fdae61", "#de2d26"),
           lwd = 2,
           bty = "n")
    dev.off()

    saveRDS(cox_model, "outputs/Clinical/nomogram_model.rds")

    model_summary <- tibble(
      C_index = cox_model$stats["C"],
      n_samples = nrow(nomogram_data),
      n_events = sum(nomogram_data$OS_status),
      genes_included = paste(top_prognostic[1:min(5, length(top_prognostic))], collapse = ", ")
    )

    export(model_summary, "outputs/Clinical/nomogram_summary.csv")

    message("Clinical nomogram completed successfully!")
    message("C-index: ", round(cox_model$stats["C"], 3))

  } else {
    message("Insufficient samples with complete clinical data for nomogram.")
    message("Available samples: ", nrow(nomogram_data))
  }

} else {
  message("Insufficient prognostic genes (need >= 3) for nomogram development.")
  message("Prognostic genes available: ", length(top_prognostic))
}
