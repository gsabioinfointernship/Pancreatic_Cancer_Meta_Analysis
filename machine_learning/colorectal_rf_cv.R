# ============================================================================
# Colorectal RF k-fold Cross-Validation (AUC correction)
# Replaces single 80/20 split with k=10 CV for publication-credible AUC
# Output: overwrites outputs/colorectal/machine_learning/colorectal_rf_metrics.csv
#         and colorectal_ml_summary.csv with CV-AUC
# ============================================================================

source("machine_learning/00_ml_functions.R")

cancer    <- "colorectal"
CV_K      <- 10   # k-fold CV folds
set.seed(42)

out_csv <- paste0("outputs/", cancer, "/machine_learning/")
out_fig <- paste0("outputs/", cancer, "/figures/machine_learning/")
dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

message("=== COLORECTAL RF k-fold CV (k=", CV_K, ") ===")

# ── 1. Load meta + hub genes
meta    <- load_meta(cancer)
hub_df  <- load_hub_genes(cancer)
candidates <- select_ml_genes(meta, hub_df, n = HUB_N)
message("ML candidates: ", nrow(candidates))

# ── 2. TCGA data
projects <- TCGA_CONFIG[[cancer]]
message("TCGA Projects: ", paste(projects, collapse = " + "))
tcga_data <- combine_tcga_ml(projects)
if (is.null(tcga_data)) stop("No TCGA data for ", cancer)

# ── 3. Prepare RF matrix (all samples)
rf_data <- prepare_rf_data(
  expr_mat    = tcga_data$expr,
  sample_type = tcga_data$sample_type,
  genes       = candidates$Symbol
)
if (is.null(rf_data)) stop("prepare_rf_data returned NULL")

X <- rf_data$X
y <- rf_data$y

message("Samples: ", nrow(X), "  Features: ", ncol(X))
message("Tumor: ", sum(y == "Tumor"), "  Normal: ", sum(y == "Normal"))

# ── 4. k-fold stratified CV
folds <- caret::createFolds(y, k = CV_K, list = TRUE, returnTrain = FALSE)

cv_aucs   <- numeric(CV_K)
cv_accs   <- numeric(CV_K)
cv_sens   <- numeric(CV_K)
cv_specs  <- numeric(CV_K)

for (i in seq_along(folds)) {
  test_idx  <- folds[[i]]
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)

  X_train <- X[train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  X_test  <- X[test_idx,  , drop = FALSE]
  y_test  <- y[test_idx]

  # Class weights to handle imbalance
  n_tumor  <- sum(y_train == "Tumor")
  n_normal <- sum(y_train == "Normal")

  rf_fit <- tryCatch(
    randomForest(
      x       = X_train,
      y       = y_train,
      ntree   = RF_TREES,
      importance = FALSE,
      classwt = c(Normal = 1, Tumor = n_normal / n_tumor)
    ),
    error = function(e) { message("  Fold ", i, " RF error: ", e$message); NULL }
  )
  if (is.null(rf_fit)) next

  pred_prob  <- predict(rf_fit, X_test, type = "prob")[, "Tumor"]
  pred_class <- predict(rf_fit, X_test, type = "class")

  roc_i <- tryCatch(
    roc(y_test, pred_prob, levels = c("Normal", "Tumor"),
        direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  cv_aucs[i]  <- if (!is.null(roc_i)) as.numeric(auc(roc_i)) else NA_real_

  cm_i <- tryCatch(
    confusionMatrix(pred_class, y_test, positive = "Tumor"),
    error = function(e) NULL
  )
  if (!is.null(cm_i)) {
    cv_accs[i]  <- cm_i$overall["Accuracy"]
    cv_sens[i]  <- cm_i$byClass["Sensitivity"]
    cv_specs[i] <- cm_i$byClass["Specificity"]
  }

  message(sprintf("  Fold %02d: AUC = %.4f  Acc = %.4f", i, cv_aucs[i], cv_accs[i]))
}

cv_auc_mean  <- mean(cv_aucs,  na.rm = TRUE)
cv_auc_sd    <- sd(cv_aucs,    na.rm = TRUE)
cv_acc_mean  <- mean(cv_accs,  na.rm = TRUE)
cv_sens_mean <- mean(cv_sens,  na.rm = TRUE)
cv_spec_mean <- mean(cv_specs, na.rm = TRUE)

message(sprintf("\n10-fold CV AUC: %.4f ± %.4f", cv_auc_mean, cv_auc_sd))
message(sprintf("10-fold CV Accuracy: %.4f", cv_acc_mean))

# ── 5. Fit final RF on ALL data (for importance export)
message("\nFitting final RF on full data for importance...")
final_rf <- randomForest(
  x          = X,
  y          = y,
  ntree      = RF_TREES,
  importance = TRUE,
  classwt    = c(Normal = 1, Tumor = sum(y == "Normal") / sum(y == "Tumor"))
)

imp <- importance(final_rf) |>
  as.data.frame() |>
  rownames_to_column("Symbol") |>
  arrange(desc(MeanDecreaseGini)) |>
  as_tibble()

# ── 6. Per-fold summary CSV
fold_summary <- data.frame(
  Fold     = seq_along(folds),
  AUC      = round(cv_aucs,  4),
  Accuracy = round(cv_accs,  4),
  Sensitivity = round(cv_sens,  4),
  Specificity = round(cv_specs, 4)
)
safe_export(fold_summary, file.path(out_csv, paste0(cancer, "_rf_cv_fold_results.csv")))
message("Saved: rf_cv_fold_results.csv")

# ── 7. Update rf_metrics.csv with CV metrics
rf_metrics_cv <- data.frame(
  Metric = c("CV_AUC_mean", "CV_AUC_sd", "CV_Accuracy", "CV_Sensitivity", "CV_Specificity",
             "n_folds", "note"),
  Value  = c(round(cv_auc_mean, 4), round(cv_auc_sd, 4),
             round(cv_acc_mean, 4), round(cv_sens_mean, 4), round(cv_spec_mean, 4),
             CV_K, "k-fold cross-validation (replaces single 80/20 holdout)")
)
safe_export(rf_metrics_cv, file.path(out_csv, paste0(cancer, "_rf_cv_metrics.csv")))
message("Saved: rf_cv_metrics.csv")

# ── 8. Update ml_summary with CV-AUC
old_summary <- read.csv(file.path(out_csv, paste0(cancer, "_ml_summary.csv")),
                        stringsAsFactors = FALSE)
old_summary$rf_auc          <- round(cv_auc_mean, 3)
old_summary$rf_auc_sd       <- round(cv_auc_sd, 3)
old_summary$rf_accuracy     <- round(cv_acc_mean, 3)
old_summary$rf_sensitivity  <- round(cv_sens_mean, 3)
old_summary$rf_specificity  <- round(cv_spec_mean, 3)
old_summary$rf_method       <- "10-fold_CV"
safe_export(old_summary, file.path(out_csv, paste0(cancer, "_ml_summary.csv")))
message("Updated: ml_summary.csv with CV-AUC")

# ── 9. Save importance plot
save_importance_plot(
  list(importance = imp),
  cancer,
  file.path(out_fig, paste0(cancer, "_rf_importance_barplot.png"))
)

# ── 10. ROC curve: use final RF on OOB predictions
oob_votes     <- final_rf$votes[, "Tumor"]
roc_oob <- tryCatch(
  roc(y, oob_votes, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE),
  error = function(e) NULL
)
if (!is.null(roc_oob)) {
  auc_oob <- as.numeric(auc(roc_oob))
  message(sprintf("OOB AUC (final model): %.4f", auc_oob))

  roc_df <- data.frame(
    specificity = roc_oob$specificities,
    sensitivity = roc_oob$sensitivities
  )
  p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = "#de2d26", linewidth = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40") +
    annotate("label", x = 0.75, y = 0.25,
             label = sprintf("CV-AUC = %.3f ± %.3f\nOOB-AUC = %.3f",
                             cv_auc_mean, cv_auc_sd, auc_oob),
             size = 4.5, fontface = "bold", fill = "white") +
    labs(
      title    = "Colorectal — Diagnostic Performance",
      subtitle = "Random Forest Classifier (10-fold CV + OOB)",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      aspect.ratio  = 1
    )
  ggsave(file.path(out_fig, paste0(cancer, "_rf_roc_cv.png")),
         plot = p_roc, width = 6, height = 6, dpi = 600, bg = "white")
  message("Saved: rf_roc_cv.png")
}

message("\n=== COLORECTAL RF CV COMPLETE ===")
message(sprintf("Final CV-AUC: %.4f ± %.4f", cv_auc_mean, cv_auc_sd))
