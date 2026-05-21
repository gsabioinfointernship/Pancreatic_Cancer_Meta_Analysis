source("validation/00_validation_functions.R")

message("=== PAN-CANCER VALIDATION ===")

out_csv <- "outputs/pancancer/validation/"
out_fig <- "outputs/pancancer/figures/validation/"
dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. Load per-cancer validation results (from 01_per_cancer_validation.R)
# ============================================================================

val_per_cancer <- lapply(setNames(CANCERS, CANCERS), function(cancer) {
  path <- paste0("outputs/", cancer, "/validation/", cancer, "_validation_results.csv")
  if (!file.exists(path)) { message("  Missing: ", path); return(NULL) }
  read.csv(path, stringsAsFactors = FALSE) |> mutate(cancer = cancer)
})
val_per_cancer <- Filter(Negate(is.null), val_per_cancer)

metrics_per_cancer <- lapply(setNames(CANCERS, CANCERS), function(cancer) {
  path <- paste0("outputs/", cancer, "/validation/", cancer, "_validation_metrics.csv")
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
})
metrics_per_cancer <- Filter(Negate(is.null), metrics_per_cancer)

if (length(val_per_cancer) == 0) {
  stop("No per-cancer validation results found. Run 01_per_cancer_validation.R first.")
}

# ============================================================================
# 2. Cross-cancer validation summary
# ============================================================================

all_val     <- bind_rows(val_per_cancer)
all_metrics <- bind_rows(metrics_per_cancer)

safe_export(all_metrics, file.path(out_csv, "pancancer_validation_metrics.csv"))
message("Validation metrics summary:")
print(all_metrics)

# ── Summary barplot: validation % per cancer
p_summary <- all_metrics |>
  pivot_longer(cols = c(concordance_pct, validation_pct),
               names_to = "Metric", values_to = "Percent") |>
  mutate(Metric = recode(Metric,
    concordance_pct = "Direction Concordance",
    validation_pct  = "Fully Validated"
  )) |>
  ggplot(aes(x = cancer, y = Percent, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            position = position_dodge(0.9), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Direction Concordance" = "#2171b5",
                                "Fully Validated"       = "#de2d26")) +
  labs(title = "Cross-Cancer Validation Rate (recount3 TCGA)",
       x = NULL, y = "Percentage (%)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank())
ggsave(file.path(out_fig, "pancancer_validation_summary.png"),
       plot = p_summary, width = 9, height = 6, dpi = 300)

# ── Spearman r barplot
p_cor <- all_metrics |>
  ggplot(aes(x = reorder(cancer, spearman_r), y = spearman_r, fill = spearman_r)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(spearman_r, 3)),
            hjust = -0.2, size = 3.5) +
  scale_fill_gradient2(low = "#2171b5", mid = "white", high = "#de2d26",
                       midpoint = 0, guide = "none") +
  coord_flip() +
  labs(title = "Meta-Analysis vs recount3 TCGA\nSpearman Correlation",
       x = NULL, y = "Spearman r") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(out_fig, "pancancer_spearman_r.png"),
       plot = p_cor, width = 7, height = 5, dpi = 300)

# ============================================================================
# 3. Consistently validated genes across cancers
# ============================================================================

message("Finding consistently validated genes...")

validated_per_cancer <- lapply(val_per_cancer, function(df) {
  df |> filter(Validated) |> pull(Symbol)
})

all_validated_symbols <- unique(unlist(validated_per_cancer))

presence_df <- tibble(Symbol = all_validated_symbols)
for (cancer in names(validated_per_cancer)) {
  presence_df[[cancer]] <- as.integer(all_validated_symbols %in% validated_per_cancer[[cancer]])
}
presence_df <- presence_df |>
  mutate(n_cancers = rowSums(across(all_of(names(validated_per_cancer))))) |>
  arrange(desc(n_cancers))

safe_export(presence_df,
            file.path(out_csv, "pancancer_validated_gene_presence.csv"))

consistently_validated <- presence_df |>
  filter(n_cancers >= MIN_CANCERS_PC)
message("Consistently validated (>=", MIN_CANCERS_PC, " cancers): ",
        nrow(consistently_validated))

safe_export(consistently_validated,
            file.path(out_csv, "pancancer_consistently_validated_genes.csv"))

# ── Bubble chart: top consistently validated genes
if (nrow(consistently_validated) > 0) {
  # Get LFC info per cancer for top genes
  top_cv <- head(consistently_validated, 40)
  bubble_df <- all_val |>
    filter(Symbol %in% top_cv$Symbol, Validated) |>
    mutate(Direction_meta = ifelse(pooled_log2FC > 0, "Up", "Down"))

  p_bubble <- ggplot(bubble_df,
                     aes(x = cancer, y = reorder(Symbol, pooled_log2FC),
                         color = Direction_meta, size = abs(pooled_log2FC))) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = c(Up = "#de2d26", Down = "#2171b5")) +
    scale_size_continuous(range = c(2, 8), name = "|Log2FC|") +
    labs(title = "Consistently Validated Genes Across Cancers",
         x = "Cancer", y = "Gene Symbol") +
    theme_bw(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 35, hjust = 1),
          plot.title   = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(out_fig, "pancancer_consistently_validated_bubble.png"),
         plot = p_bubble,
         width  = 10,
         height = min(48, max(6, nrow(top_cv) * 0.35 + 3)),
         dpi    = 300)
}

# ── Cross-cancer LFC correlation heatmap (spearman r matrix)
message("Building cross-cancer LFC correlation matrix...")

lfc_wide <- all_val |>
  select(Symbol, cancer, pooled_log2FC) |>
  pivot_wider(names_from = cancer, values_from = pooled_log2FC)

if (ncol(lfc_wide) > 2) {
  lfc_mat <- lfc_wide |> select(-Symbol) |> as.matrix()
  rownames(lfc_mat) <- lfc_wide$Symbol
  lfc_mat <- lfc_mat[complete.cases(lfc_mat), ]

  if (nrow(lfc_mat) > 10) {
    cor_mat <- cor(lfc_mat, method = "spearman", use = "pairwise.complete.obs")

    png(file.path(out_fig, "pancancer_lfc_correlation_matrix.png"),
        width = 7, height = 6, units = "in", res = 300)
    pheatmap(
      cor_mat,
      display_numbers = TRUE,
      number_format   = "%.2f",
      color           = colorRampPalette(c("#2171b5", "white", "#de2d26"))(100),
      main            = "Cross-Cancer LFC Correlation (Spearman)",
      fontsize         = 11,
      fontsize_number  = 9
    )
    dev.off()
  }
}

message("\nPan-cancer validation complete.")
