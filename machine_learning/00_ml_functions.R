library(glmnet)
library(survival)
library(survminer)
library(randomForest)
library(pROC)
library(caret)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# ── Cancer → TCGA project mapping
TCGA_CONFIG <- list(
  colorectal = c("COAD", "READ"),
  esophagus  = c("ESCA"),
  kidney     = c("KIRC", "KIRP"),
  liver      = c("LIHC"),
  pancreatic = c("PAAD")
)

CANCERS     <- names(TCGA_CONFIG)
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1
HUB_N       <- 100      # top hub genes as feature candidates
LASSO_ALPHA <- 0.9       # 0.9 = Elastic Net (more stable than pure LASSO)
CV_FOLDS    <- 10        # cross-validation folds for glmnet
MIN_EVENTS  <- 15        # minimum deaths for LASSO-Cox
MIN_NORMAL  <- 5         # minimum normal samples for RF classifier
RF_TREES    <- 500       # random forest trees
TRAIN_FRAC  <- 0.8       # train/test split fraction
GLMNET_MAXIT <- 1e6      # higher iterations for convergence
GLMNET_THRESH <- 1e-6    # slightly relaxed threshold
set.seed(42)

# ── Load meta-analysis results
load_meta <- function(cancer) {
  path <- paste0("outputs/", cancer, "/meta_analysis/", cancer, "_meta_analysis_results.rds")
  readRDS(path)
}

# ── Load hub genes
load_hub_genes <- function(cancer) {
  path <- paste0("outputs/", cancer, "/network/", cancer, "_hub_genes.csv")
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}

# ── Download TCGA (all samples: tumor + normal) via curatedTCGAData
download_tcga_ml <- function(project_id) {
  message("    Downloading: ", project_id)
  tryCatch(
    curatedTCGAData(
      diseaseCode = project_id,
      assays      = "RNASeq2GeneNorm",
      version     = "2.1.1",
      dry.run     = FALSE
    ),
    error = function(e) { message("    Error: ", e$message); NULL }
  )
}

# ── Extract expression + labels + OS clinical from MAE
extract_ml_data <- function(mae, project_id) {
  assay_name <- grep("RNASeq2GeneNorm", names(mae), value = TRUE)[1]
  if (is.na(assay_name)) return(NULL)

  expr_mat  <- assay(mae[[assay_name]])
  type_code <- as.integer(substr(colnames(expr_mat), 14, 15))
  sample_type <- ifelse(type_code < 10, "Tumor", "Normal")

  # log2 transform if needed
  if (max(expr_mat, na.rm = TRUE) > 100) expr_mat <- log2(expr_mat + 1)

  # Remove non-finite rows
  expr_mat <- expr_mat[apply(expr_mat, 1, function(x) all(is.finite(x))), ]

  # Clinical
  clin <- tryCatch(
    as.data.frame(colData(mae), stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )

  pid_col      <- grep("patientID|submitter_id", colnames(clin), ignore.case = TRUE, value = TRUE)[1]
  vital_col    <- grep("vital_status",           colnames(clin), ignore.case = TRUE, value = TRUE)[1]
  death_col    <- grep("days_to_death",           colnames(clin), ignore.case = TRUE, value = TRUE)[1]
  followup_col <- grep("days_to_last_follow",     colnames(clin), ignore.case = TRUE, value = TRUE)[1]
  stage_col    <- grep("pathologic_stage|tumor_stage|ajcc_pathologic_tumor_stage",
                       colnames(clin), ignore.case = TRUE, value = TRUE)[1]
  age_col      <- grep("age_at_initial|age_at_diagnosis",
                       colnames(clin), ignore.case = TRUE, value = TRUE)[1]

  os_df <- NULL
  if (!is.na(pid_col) && !is.na(vital_col) && nrow(clin) > 0) {
    os_df <- clin |>
      as_tibble() |>
      mutate(patient_id = .data[[pid_col]]) |>
      select(patient_id,
             any_of(c(vital_col, death_col, followup_col, stage_col, age_col))) |>
      filter(!duplicated(patient_id)) |>
      mutate(
        OS_status = as.integer(tolower(.data[[vital_col]]) %in% c("dead", "deceased", "1")),
        OS_days   = {
          d <- if (!is.na(death_col)    && death_col    %in% colnames(pick(everything()))) as.numeric(.data[[death_col]])    else NA_real_
          f <- if (!is.na(followup_col) && followup_col %in% colnames(pick(everything()))) as.numeric(.data[[followup_col]]) else NA_real_
          dplyr::coalesce(d, f)
        }
      ) |>
      filter(!is.na(OS_days), is.finite(OS_days), OS_days > 0)
  }

  list(
    expr        = expr_mat,
    sample_type = sample_type,
    clinical    = os_df,
    project     = project_id
  )
}

# ── Combine multiple TCGA projects
combine_tcga_ml <- function(project_ids) {
  all_data <- lapply(project_ids, function(p) {
    mae <- download_tcga_ml(p)
    if (is.null(mae)) return(NULL)
    extract_ml_data(mae, p)
  })
  all_data <- Filter(Negate(is.null), all_data)
  if (length(all_data) == 0) return(NULL)

  common_genes <- Reduce(intersect, lapply(all_data, function(x) rownames(x$expr)))
  expr         <- do.call(cbind, lapply(all_data, function(x) x$expr[common_genes, ]))
  sample_type  <- unlist(lapply(all_data, function(x) x$sample_type))
  clinical     <- bind_rows(lapply(all_data, function(x) x$clinical)) |>
    filter(!is.na(patient_id), !duplicated(patient_id))

  list(expr = expr, sample_type = sample_type, clinical = clinical)
}

# ── Select candidate genes for ML
select_ml_genes <- function(meta_df, hub_df, n = HUB_N) {
  sig <- meta_df |>
    filter(rem_padj < PADJ_CUTOFF, abs(pooled_log2FC) >= LFC_CUTOFF) |>
    select(Symbol, pooled_log2FC, rem_padj)

  if (!is.null(hub_df) && nrow(hub_df) > 0) {
    hub_df |> inner_join(sig, by = "Symbol") |> arrange(desc(Degree)) |> head(n)
  } else {
    sig |> arrange(rem_padj) |> head(n)
  }
}

# ============================================================================
# LASSO-COX RISK SCORE
# ============================================================================

# ── Prepare feature matrix for LASSO-Cox (tumor samples with OS data)
prepare_lasso_data <- function(expr_mat, sample_type, clinical, genes) {
  # Tumor samples only
  tumor_idx   <- which(sample_type == "Tumor")
  tumor_expr  <- expr_mat[, tumor_idx, drop = FALSE]
  patient_ids <- substr(colnames(tumor_expr), 1, 12)

  # Join with OS
  df <- tibble(
    sample_id  = colnames(tumor_expr),
    patient_id = patient_ids
  ) |>
    left_join(clinical |> select(patient_id, OS_days, OS_status),
              by = "patient_id") |>
    filter(!is.na(OS_days), !is.na(OS_status), OS_days > 0) |>
    filter(!duplicated(patient_id))

  if (nrow(df) < 20) {
    message("    Too few patients with OS data (n=", nrow(df), ")")
    return(NULL)
  }
  if (sum(df$OS_status) < MIN_EVENTS) {
    message("    Too few death events (n=", sum(df$OS_status), " < ", MIN_EVENTS, ")")
    return(NULL)
  }

  # Feature matrix: patients × genes
  avail_genes <- intersect(genes, rownames(tumor_expr))
  X <- t(tumor_expr[avail_genes, df$sample_id, drop = FALSE])

  # Scale features
  X <- scale(X)
  X[!is.finite(X)] <- 0

  y <- Surv(df$OS_days, df$OS_status)

  list(X = X, y = y, df = df, genes = avail_genes)
}

# ── Run LASSO-Cox with cross-validation
run_lasso_cox <- function(lasso_data) {
  X <- lasso_data$X
  y <- lasso_data$y

  message("    Running LASSO-Cox: ", nrow(X), " patients, ",
          ncol(X), " genes, ", sum(lasso_data$df$OS_status), " events")

  # 10-fold CV
  cv_fit <- tryCatch(
    cv.glmnet(X, y, family = "cox", alpha = LASSO_ALPHA,
              nfolds = CV_FOLDS, type.measure = "C",
              maxit = GLMNET_MAXIT, thresh = GLMNET_THRESH),
    error = function(e) { message("    glmnet error: ", e$message); NULL }
  )
  if (is.null(cv_fit)) return(NULL)

  # Use lambda.min (best CV performance)
  lambda_best <- cv_fit$lambda.min

  # Extract non-zero coefficients
  coefs <- as.matrix(coef(cv_fit, s = lambda_best))
  selected_coefs <- coefs[coefs[, 1] != 0, , drop = FALSE]

  n_selected <- nrow(selected_coefs)
  message("    LASSO selected genes: ", n_selected,
          "  (lambda.min = ", round(lambda_best, 4), ")")

  if (n_selected == 0) {
    # Try lambda.1se
    lambda_best <- cv_fit$lambda.1se
    coefs       <- as.matrix(coef(cv_fit, s = lambda_best))
    selected_coefs <- coefs[coefs[, 1] != 0, , drop = FALSE]
    n_selected  <- nrow(selected_coefs)
    message("    Fallback lambda.1se — selected: ", n_selected)
  }

  if (n_selected == 0) {
    message("    No genes selected by LASSO")
    return(NULL)
  }

  list(
    cv_fit         = cv_fit,
    lambda         = lambda_best,
    selected_genes = rownames(selected_coefs),
    coefs          = selected_coefs,
    n_selected     = n_selected
  )
}

# ── Compute risk score and evaluate
compute_and_evaluate_risk <- function(lasso_result, lasso_data, cancer) {
  X  <- lasso_data$X
  df <- lasso_data$df

  # Risk score: linear combination of selected gene expressions × coefficients
  sel_genes <- lasso_result$selected_genes
  X_sel     <- X[, sel_genes, drop = FALSE]
  risk_score <- as.numeric(X_sel %*% lasso_result$coefs[sel_genes, 1])

  df$risk_score <- risk_score
  df$risk_group <- factor(
    ifelse(risk_score >= median(risk_score), "High", "Low"),
    levels = c("Low", "High")
  )

  # C-index
  cx_full <- tryCatch(
    coxph(Surv(OS_days, OS_status) ~ risk_score, data = df),
    error = function(e) NULL
  )
  c_index <- if (!is.null(cx_full)) summary(cx_full)$concordance[1] else NA_real_
  message("    C-index: ", round(c_index, 3))

  # Log-rank test
  lrt <- tryCatch(
    survdiff(Surv(OS_days, OS_status) ~ risk_group, data = df),
    error = function(e) NULL
  )
  km_pval <- if (!is.null(lrt)) 1 - pchisq(lrt$chisq, df = 1) else NA_real_

  list(df = df, c_index = c_index, km_pval = km_pval, cx_model = cx_full)
}

# ============================================================================
# RANDOM FOREST CLASSIFIER (Tumor vs Normal)
# ============================================================================

# ── Prepare RF data: all samples (tumor + normal), feature = hub gene expression
prepare_rf_data <- function(expr_mat, sample_type, genes) {
  n_tumor  <- sum(sample_type == "Tumor")
  n_normal <- sum(sample_type == "Normal")
  message("    Tumor: ", n_tumor, "  Normal: ", n_normal)

  if (n_normal < MIN_NORMAL) {
    message("    Too few normal samples (n=", n_normal, ") — skipping RF")
    return(NULL)
  }

  avail_genes <- intersect(genes, rownames(expr_mat))
  if (length(avail_genes) < 5) {
    message("    Too few genes in expression matrix")
    return(NULL)
  }

  X <- t(expr_mat[avail_genes, , drop = FALSE])
  y <- factor(sample_type, levels = c("Normal", "Tumor"))

  # Remove columns with zero variance
  vars <- apply(X, 2, var, na.rm = TRUE)
  X    <- X[, vars > 0, drop = FALSE]
  X[!is.finite(X)] <- 0

  message("    RF features: ", ncol(X), " genes, ", nrow(X), " samples")
  list(X = X, y = y, genes = colnames(X))
}

# ── Train RF with train/test split + CV
run_rf_classifier <- function(rf_data) {
  X <- rf_data$X
  y <- rf_data$y

  # Stratified train/test split
  set.seed(42)
  idx_tumor  <- which(y == "Tumor")
  idx_normal <- which(y == "Normal")
  train_idx  <- c(
    sample(idx_tumor,  floor(length(idx_tumor)  * TRAIN_FRAC)),
    sample(idx_normal, floor(length(idx_normal) * TRAIN_FRAC))
  )
  test_idx   <- setdiff(seq_len(nrow(X)), train_idx)

  X_train <- X[train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  X_test  <- X[test_idx,  , drop = FALSE]
  y_test  <- y[test_idx]

  message("    Training RF: ", length(train_idx), " train, ", length(test_idx), " test")

  rf_fit <- tryCatch(
    randomForest(
      x          = X_train,
      y          = y_train,
      ntree      = RF_TREES,
      importance = TRUE,
      classwt    = c(Normal = 1, Tumor = length(idx_normal) / length(idx_tumor))
    ),
    error = function(e) { message("    RF error: ", e$message); NULL }
  )
  if (is.null(rf_fit)) return(NULL)

  # Predict test set
  pred_prob  <- predict(rf_fit, X_test, type = "prob")[, "Tumor"]
  pred_class <- predict(rf_fit, X_test, type = "class")

  # ROC + AUC
  roc_obj <- tryCatch(
    roc(y_test, pred_prob, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  auc_val <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA_real_
  message("    RF AUC: ", round(auc_val, 3))

  # Confusion matrix
  cm <- confusionMatrix(pred_class, y_test, positive = "Tumor")

  # Feature importance (top genes by MeanDecreaseGini)
  imp <- importance(rf_fit) |>
    as.data.frame() |>
    rownames_to_column("Symbol") |>
    arrange(desc(MeanDecreaseGini)) |>
    as_tibble()

  list(
    rf_fit    = rf_fit,
    roc       = roc_obj,
    auc       = auc_val,
    cm        = cm,
    accuracy  = cm$overall["Accuracy"],
    sens      = cm$byClass["Sensitivity"],
    spec      = cm$byClass["Specificity"],
    importance = imp,
    y_test    = y_test,
    pred_prob = pred_prob
  )
}

# ============================================================================
# FIGURES
# ============================================================================

# ── LASSO CV curve
save_lasso_cv_plot <- function(cv_fit, cancer, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  # Use higher DPI for publication
  png(path, width = 6, height = 5, units = "in", res = 600)

  # Adjust margins: bottom, left, top, right
  par(mar = c(5, 5, 5, 2), mgp = c(3, 1, 0), las = 1, font.main = 2)

  # Plot CV curve
  plot(cv_fit)

  # Add custom title that doesn't overlap with non-zero counts
  title(main = paste(tools::toTitleCase(cancer), "— LASSO-Cox CV"),
        line = 3, cex.main = 1.1)

  # Highlight lambda.min
  abline(v = log(cv_fit$lambda.min), col = "#de2d26", lty = 2, lwd = 1.5)

  # Better legend
  legend("topright",
         legend = paste0("log(lambda.min) = ", round(log(cv_fit$lambda.min), 3)),
         col = "#de2d26", lty = 2, lwd = 1.5, cex = 0.8, bty = "n")

  dev.off()
}

# ── Risk score KM plot
save_risk_km_plot <- function(risk_df, cancer, c_index, km_pval, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (is.null(risk_df)) return(invisible(NULL))

  fit <- tryCatch(
    survfit(Surv(OS_days, OS_status) ~ risk_group, data = risk_df),
    error = function(e) NULL
  )
  if (is.null(fit)) return(invisible(NULL))

  # Cleaner subtitle
  subtitle_txt <- sprintf("Log-rank p = %.3g  |  C-index = %.3f",
                           ifelse(is.na(km_pval), 1, km_pval),
                           ifelse(is.na(c_index), 0, c_index))

  p <- ggsurvplot(
    fit, data = risk_df,
    palette      = c("#3182bd", "#de2d26"), # Better blue/red
    title        = paste(tools::toTitleCase(cancer), "— Risk Model Survival"),
    subtitle     = subtitle_txt,
    xlab         = "Time (Days)",
    ylab         = "Overall Survival Probability",
    legend.title = "Risk Group",
    legend.labs  = c("Low Risk", "High Risk"),
    risk.table   = TRUE,
    risk.table.y.text = FALSE,
    risk.table.height = 0.25,
    ggtheme      = theme_classic(base_size = 12) +
                   theme(plot.title = element_text(face = "bold", hjust = 0.5),
                         plot.subtitle = element_text(hjust = 0.5)),
    conf.int     = TRUE,
    conf.int.alpha = 0.1,
    censor.size  = 3
  )

  tryCatch({
    png(path, width = 7, height = 7, units = "in", res = 600)
    print(p)
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("    KM save error: ", e$message)
  })
  invisible(NULL)
}

# ── Risk score distribution (density plot)
save_risk_distribution <- function(risk_df, cancer, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (is.null(risk_df)) return(invisible(NULL))

  p <- ggplot(risk_df, aes(x = risk_score, fill = risk_group)) +
    geom_histogram(aes(y = after_stat(density)), bins = 35, alpha = 0.5, position = "identity") +
    geom_density(aes(color = risk_group), linewidth = 1, fill = NA) +
    geom_vline(xintercept = median(risk_df$risk_score),
               linetype = "dashed", color = "gray20", linewidth = 0.8) +
    scale_fill_manual(values  = c("Low" = "#3182bd", "High" = "#de2d26")) +
    scale_color_manual(values = c("Low" = "#3182bd", "High" = "#de2d26")) +
    labs(
      title    = paste(tools::toTitleCase(cancer), "— Risk Score Distribution"),
      subtitle = "Median-split used for group definition",
      x        = "LASSO-Cox Risk Score",
      y        = "Density",
      fill     = "Risk Group",
      color    = "Risk Group"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "right"
    )

  ggsave(path, plot = p, width = 7, height = 5, dpi = 600, bg = "white")
  invisible(p)
}

# ── LASSO selected gene coefficients barplot
save_coef_plot <- function(coefs_df, cancer, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (is.null(coefs_df) || nrow(coefs_df) == 0) return(invisible(NULL))

  p <- ggplot(coefs_df, aes(x = reorder(Symbol, Coefficient),
                             y = Coefficient,
                             fill = Coefficient > 0)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#de2d26", "FALSE" = "#3182bd"),
                      labels = c("TRUE" = "High risk", "FALSE" = "Protective"),
                      name   = "Role") +
    labs(
      title = paste(tools::toTitleCase(cancer), "— Prognostic Features"),
      subtitle = "LASSO-Cox Non-zero Coefficients",
      x     = NULL,
      y     = "Coefficient"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.text.y = element_text(face = "italic")
    )

  ht <- max(4, nrow(coefs_df) * 0.3 + 2)
  ggsave(path, plot = p, width = 6, height = min(15, ht), dpi = 600, bg = "white")
  invisible(p)
}

# ── RF ROC curve
save_roc_plot <- function(rf_result, cancer, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (is.null(rf_result) || is.null(rf_result$roc)) return(invisible(NULL))

  roc_df <- data.frame(
    specificity = rf_result$roc$specificities,
    sensitivity = rf_result$roc$sensitivities
  )

  p <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = "#de2d26", linewidth = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray40") +
    annotate("label", x = 0.75, y = 0.2,
             label = sprintf("AUC = %.3f", rf_result$auc),
             size = 5, fontface = "bold", fill = "white") +
    labs(
      title = paste(tools::toTitleCase(cancer), "— Diagnostic Performance"),
      subtitle = "Random Forest Classifier (Tumor vs Normal)",
      x     = "False Positive Rate (1 - Specificity)",
      y     = "True Positive Rate (Sensitivity)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      aspect.ratio = 1
    )

  ggsave(path, plot = p, width = 6, height = 6, dpi = 600, bg = "white")
  invisible(p)
}

# ── RF feature importance barplot (top 20)
save_importance_plot <- function(rf_result, cancer, path, n = 20) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (is.null(rf_result) || nrow(rf_result$importance) == 0) return(invisible(NULL))

  plot_df <- rf_result$importance |>
    head(n) |>
    mutate(Symbol = factor(Symbol, levels = rev(Symbol)))

  p <- ggplot(plot_df, aes(x = MeanDecreaseGini, y = Symbol, fill = MeanDecreaseGini)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_gradient(low = "#deebf7", high = "#084594", guide = "none") +
    labs(
      title = paste(tools::toTitleCase(cancer), "— Diagnostic Hub Genes"),
      subtitle = paste("Top", n, "Predictors by Mean Decrease Gini"),
      x     = "Feature Importance (Mean Decrease Gini)",
      y     = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.text.y = element_text(face = "italic")
    )

  ht <- max(4, n * 0.3 + 2)
  ggsave(path, plot = p, width = 7, height = min(12, ht), dpi = 600, bg = "white")
  invisible(p)
}

# ── Safe CSV export
safe_export <- function(df, path) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- sapply(df[[nm]], function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        tryCatch(paste(unique(rapply(x, as.character, how = "unlist")), collapse = "; "),
                 error = function(e) NA_character_)
      }, USE.NAMES = FALSE)
    }
  }
  write.csv(df, path, row.names = FALSE)
}
