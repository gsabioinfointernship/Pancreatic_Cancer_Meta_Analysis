library(caret)
library(randomForest)
library(e1071)
library(glmnet)
library(xgboost)
library(pROC)
library(tidyverse)
library(rio)

dir.create("outputs/MachineLearning", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/MachineLearning", showWarnings = FALSE, recursive = TRUE)

sig_genes <- import("outputs/MetaVolcanoR/significant_genes.csv")

feature_sets <- list(
  top50 = sig_genes %>% arrange(randomP) %>% head(50) %>% pull(Gene_ID),
  top100 = sig_genes %>% arrange(randomP) %>% head(100) %>% pull(Gene_ID),
  top200 = sig_genes %>% arrange(randomP) %>% head(200) %>% pull(Gene_ID)
)

geo_data <- list()
geo_ids <- c("GSE130688", "GSE136569", "GSE171485", "GSE196009",
             "GSE211398", "GSE280271", "GSE293744")

for(geo_id in geo_ids) {
  metadata <- import(paste0("data/", geo_id, "_metadata.csv"))
  counts <- import(paste0("data/", geo_id, "_raw_counts.tsv"))

  counts_matrix <- counts %>% column_to_rownames("GeneID") %>% as.matrix()

  vst_counts <- DESeq2::vst(counts_matrix + 1)

  geo_data[[geo_id]] <- list(
    expression = vst_counts,
    labels = metadata$condition
  )
}

tcga_normalized <- readRDS("data/TCGA/TCGA_PAAD_normalized.rds")
tcga_clinical <- import("data/TCGA/TCGA_PAAD_clinical.csv")

results_all <- tibble()

for(feature_name in names(feature_sets)) {

  genes <- feature_sets[[feature_name]]

  geo_expr_list <- lapply(geo_data, function(x) {
    common_genes <- intersect(genes, rownames(x$expression))
    t(x$expression[common_genes, ])
  })

  geo_expr <- do.call(rbind, geo_expr_list)
  geo_labels <- unlist(lapply(geo_data, function(x) x$labels))
  geo_labels <- factor(ifelse(geo_labels %in% c("normal", "Normal"), "Normal", "Tumor"))

  train_data <- as.data.frame(geo_expr)
  train_data$outcome <- geo_labels

  common_genes <- intersect(genes, rownames(tcga_normalized))
  tcga_expr <- t(tcga_normalized[common_genes, tcga_clinical$sample_id])

  test_data <- as.data.frame(tcga_expr)
  test_labels <- factor(tcga_clinical$sample_type, levels = c("Normal", "Tumor"))
  test_data$outcome <- test_labels

  train_data <- train_data[complete.cases(train_data), ]
  test_data <- test_data[complete.cases(test_data), ]

  train_control <- trainControl(
    method = "cv",
    number = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )

  set.seed(123)
  rf_model <- train(
    outcome ~ .,
    data = train_data,
    method = "rf",
    trControl = train_control,
    metric = "ROC",
    ntree = 500
  )

  set.seed(123)
  svm_model <- train(
    outcome ~ .,
    data = train_data,
    method = "svmRadial",
    trControl = train_control,
    metric = "ROC",
    preProcess = c("center", "scale")
  )

  set.seed(123)
  enet_model <- train(
    outcome ~ .,
    data = train_data,
    method = "glmnet",
    trControl = train_control,
    metric = "ROC",
    family = "binomial"
  )

  xgb_train <- xgb.DMatrix(
    data = as.matrix(train_data[, -ncol(train_data)]),
    label = as.numeric(train_data$outcome) - 1
  )

  set.seed(123)
  xgb_model <- xgb.cv(
    data = xgb_train,
    nrounds = 100,
    nfold = 10,
    objective = "binary:logistic",
    eta = 0.3,
    max_depth = 6,
    early_stopping_rounds = 10,
    verbose = 0
  )

  xgb_final <- xgboost(
    data = xgb_train,
    nrounds = xgb_model$best_iteration,
    objective = "binary:logistic",
    eta = 0.3,
    max_depth = 6,
    verbose = 0
  )

  models <- list(
    RF = rf_model,
    SVM = svm_model,
    ElasticNet = enet_model,
    XGBoost = xgb_final
  )

  for(model_name in names(models)) {

    if(model_name == "XGBoost") {
      pred_train <- predict(models[[model_name]], as.matrix(train_data[, -ncol(train_data)]))
      pred_test <- predict(models[[model_name]], as.matrix(test_data[, -ncol(test_data)]))
    } else {
      pred_train <- predict(models[[model_name]], train_data, type = "prob")[, "Tumor"]
      pred_test <- predict(models[[model_name]], test_data, type = "prob")[, "Tumor"]
    }

    roc_train <- roc(train_data$outcome, pred_train, levels = c("Normal", "Tumor"), direction = "<")
    roc_test <- roc(test_data$outcome, pred_test, levels = c("Normal", "Tumor"), direction = "<")

    pred_class_train <- ifelse(pred_train > 0.5, "Tumor", "Normal")
    pred_class_test <- ifelse(pred_test > 0.5, "Tumor", "Normal")

    cm_train <- confusionMatrix(factor(pred_class_train, levels = c("Normal", "Tumor")),
                                train_data$outcome)
    cm_test <- confusionMatrix(factor(pred_class_test, levels = c("Normal", "Tumor")),
                               test_data$outcome)

    results_all <- bind_rows(results_all, tibble(
      Feature_Set = feature_name,
      Model = model_name,
      Dataset = "Training (GEO)",
      AUC = as.numeric(roc_train$auc),
      Accuracy = cm_train$overall["Accuracy"],
      Sensitivity = cm_train$byClass["Sensitivity"],
      Specificity = cm_train$byClass["Specificity"],
      PPV = cm_train$byClass["Pos Pred Value"],
      NPV = cm_train$byClass["Neg Pred Value"]
    ))

    results_all <- bind_rows(results_all, tibble(
      Feature_Set = feature_name,
      Model = model_name,
      Dataset = "Validation (TCGA)",
      AUC = as.numeric(roc_test$auc),
      Accuracy = cm_test$overall["Accuracy"],
      Sensitivity = cm_test$byClass["Sensitivity"],
      Specificity = cm_test$byClass["Specificity"],
      PPV = cm_test$byClass["Pos Pred Value"],
      NPV = cm_test$byClass["Neg Pred Value"]
    ))
  }
}

export(results_all, "outputs/MachineLearning/model_performance.csv")

best_features <- "top100"
genes <- feature_sets[[best_features]]

geo_expr_list <- lapply(geo_data, function(x) {
  common_genes <- intersect(genes, rownames(x$expression))
  t(x$expression[common_genes, ])
})

geo_expr <- do.call(rbind, geo_expr_list)
geo_labels <- factor(unlist(lapply(geo_data, function(x)
  ifelse(x$labels %in% c("normal", "Normal"), "Normal", "Tumor"))))

train_data <- as.data.frame(geo_expr)
train_data$outcome <- geo_labels
train_data <- train_data[complete.cases(train_data), ]

common_genes <- intersect(genes, rownames(tcga_normalized))
tcga_expr <- t(tcga_normalized[common_genes, tcga_clinical$sample_id])

test_data <- as.data.frame(tcga_expr)
test_labels <- factor(tcga_clinical$sample_type, levels = c("Normal", "Tumor"))
test_data$outcome <- test_labels
test_data <- test_data[complete.cases(test_data), ]

train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE)

set.seed(123)
final_rf <- train(outcome ~ ., data = train_data, method = "rf",
                  trControl = train_control, ntree = 500)

pred_train_final <- predict(final_rf, train_data, type = "prob")[, "Tumor"]
pred_test_final <- predict(final_rf, test_data, type = "prob")[, "Tumor"]

roc_train_final <- roc(train_data$outcome, pred_train_final,
                       levels = c("Normal", "Tumor"), direction = "<")
roc_test_final <- roc(test_data$outcome, pred_test_final,
                      levels = c("Normal", "Tumor"), direction = "<")

png("figures/MachineLearning/ROC_curves.png", width = 10, height = 8, units = "in", res = 300)
plot(roc_train_final, col = "#2171b5", lwd = 2, main = "ROC Curves - Random Forest")
plot(roc_test_final, col = "#de2d26", lwd = 2, add = TRUE)
legend("bottomright",
       legend = c(paste0("Training (AUC = ", round(roc_train_final$auc, 3), ")"),
                  paste0("Validation (AUC = ", round(roc_test_final$auc, 3), ")")),
       col = c("#2171b5", "#de2d26"), lwd = 2)
dev.off()

importance <- varImp(final_rf)$importance %>%
  rownames_to_column("Gene") %>%
  arrange(desc(Overall)) %>%
  head(20)

export(importance, "outputs/MachineLearning/feature_importance.csv")

p_importance <- ggplot(importance, aes(x = reorder(Gene, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "#2171b5") +
  coord_flip() +
  labs(title = "Top 20 Important Features (Random Forest)",
       x = "Gene", y = "Importance") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/MachineLearning/feature_importance.png", p_importance,
       width = 8, height = 6, dpi = 300)

message("Machine learning analysis completed!")
message("Best validation AUC: ", round(max(results_all %>%
  filter(Dataset == "Validation (TCGA)") %>% pull(AUC)), 3))
