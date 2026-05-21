source("survival_analysis/00_survival_functions.R")

message("=== PAN-CANCER SURVIVAL ANALYSIS ===")

out_csv <- "outputs/pancancer/survival/"
out_fig <- "outputs/pancancer/figures/survival/"
dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. Load per-cancer Cox results
# ============================================================================

cox_per_cancer <- lapply(setNames(CANCERS, CANCERS), function(cancer) {
  path <- paste0("outputs/", cancer, "/survival/", cancer, "_univariate_cox.csv")
  if (!file.exists(path)) { message("  Missing: ", path); return(NULL) }
  read.csv(path, stringsAsFactors = FALSE) |> mutate(cancer = cancer)
})
cox_per_cancer <- Filter(Negate(is.null), cox_per_cancer)

metrics_per_cancer <- lapply(setNames(CANCERS, CANCERS), function(cancer) {
  path <- paste0("outputs/", cancer, "/survival/", cancer, "_survival_metrics.csv")
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
})
metrics_per_cancer <- Filter(Negate(is.null), metrics_per_cancer)

if (length(cox_per_cancer) == 0) {
  stop("No per-cancer Cox results found. Run 01_per_cancer_survival.R first.")
}

message("Cancers with results: ", paste(names(cox_per_cancer), collapse = ", "))

all_cox     <- bind_rows(cox_per_cancer)
all_metrics <- bind_rows(metrics_per_cancer)

safe_export(all_metrics, file.path(out_csv, "pancancer_survival_metrics.csv"))
message("\nSurvival metrics per cancer:")
print(all_metrics |> select(cancer, n_tested, n_events, n_prognostic,
                              n_high_risk, n_protective, top_gene, top_HR))

# ============================================================================
# 2. Identify consistently prognostic genes (>= MIN_CANCERS)
# ============================================================================

message("\nFinding consistently prognostic genes (FDR < ", COX_PADJ, " in >= ",
        MIN_CANCERS, " cancers)...")

prog_per_cancer <- lapply(cox_per_cancer, function(df) {
  df |> filter(padj < COX_PADJ) |> pull(Symbol)
})

all_prog_symbols <- unique(unlist(prog_per_cancer))
message("Total unique prognostic genes: ", length(all_prog_symbols))

if (length(all_prog_symbols) > 0) {
  # Presence matrix
  presence_df <- tibble(Symbol = all_prog_symbols)
  for (cancer in names(prog_per_cancer)) {
    presence_df[[cancer]] <- as.integer(all_prog_symbols %in% prog_per_cancer[[cancer]])
  }
  presence_df <- presence_df |>
    mutate(n_cancers = rowSums(across(all_of(names(prog_per_cancer))))) |>
    arrange(desc(n_cancers))

  safe_export(presence_df, file.path(out_csv, "pancancer_prognostic_gene_presence.csv"))

  # Consistently prognostic across >= MIN_CANCERS
  consistent <- presence_df |> filter(n_cancers >= MIN_CANCERS)
  message("Consistently prognostic (>= ", MIN_CANCERS, " cancers): ", nrow(consistent))
  safe_export(consistent, file.path(out_csv, "pancancer_consistent_prognostic_genes.csv"))
} else {
  consistent  <- tibble(Symbol = character())
  presence_df <- tibble(Symbol = character())
}

# ============================================================================
# 3. HR matrix across cancers
# ============================================================================

message("\nBuilding HR matrix...")

hr_wide <- all_cox |>
  select(Symbol, cancer, HR) |>
  pivot_wider(names_from = cancer, values_from = HR)

pval_wide <- all_cox |>
  select(Symbol, cancer, padj) |>
  pivot_wider(names_from = cancer, values_from = padj,
              names_prefix = "padj_")

hr_annotated <- hr_wide |>
  left_join(pval_wide, by = "Symbol") |>
  left_join(presence_df |> select(Symbol, n_cancers), by = "Symbol")

safe_export(hr_annotated, file.path(out_csv, "pancancer_hr_matrix.csv"))

# ============================================================================
# 4. Summary barplot — n prognostic genes per cancer
# ============================================================================

p_bar <- all_metrics |>
  filter(!is.na(n_prognostic)) |>
  ggplot(aes(x = reorder(cancer, n_prognostic), y = n_prognostic, fill = cancer)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = n_prognostic), hjust = -0.2, size = 4) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Prognostic Hub Genes Per Cancer (FDR < 0.05)",
    x     = NULL,
    y     = "Number of Prognostic Genes"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(out_fig, "pancancer_prognostic_gene_counts.png"),
       plot = p_bar, width = 7, height = 5, dpi = 300)

# ============================================================================
# 5. Dot plot — prognostic genes × cancer (top consistent genes)
# ============================================================================

if (nrow(consistent) > 0) {
  top_consistent <- head(consistent, 40)

  dot_df <- all_cox |>
    filter(Symbol %in% top_consistent$Symbol, padj < COX_PADJ) |>
    mutate(
      log10p    = -log10(padj),
      Direction = ifelse(HR > 1, "High risk (HR>1)", "Protective (HR<1)")
    )

  if (nrow(dot_df) > 0) {
    p_dot <- ggplot(dot_df,
                    aes(x = cancer,
                        y = reorder(Symbol, HR),
                        color = Direction,
                        size  = log10p)) +
      geom_point(alpha = 0.85) +
      scale_color_manual(values = c(
        "High risk (HR>1)"   = "#de2d26",
        "Protective (HR<1)"  = "#2171b5"
      )) +
      scale_size_continuous(range = c(2, 8), name = "-Log10(FDR)") +
      labs(
        title  = "Consistently Prognostic Hub Genes Across Cancers",
        x      = "Cancer",
        y      = "Gene Symbol",
        color  = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x  = element_text(angle = 35, hjust = 1),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
      )

    ht <- min(20, max(5, nrow(top_consistent) * 0.3 + 3))
    ggsave(file.path(out_fig, "pancancer_consistent_prognostic_dot.png"),
           plot = p_dot, width = 9, height = ht, dpi = 300)
    message("Dot plot saved: ", nrow(dot_df), " gene-cancer pairs")
  }
}

# ============================================================================
# 6. HR heatmap across cancers (log2 scale)
# ============================================================================

message("Building HR heatmap...")

# Use all genes prognostic in >= 1 cancer, show HR for all cancers
if (length(all_prog_symbols) > 0) {
  top_genes_hm <- if (nrow(consistent) > 0) {
    # Prefer consistently prognostic genes
    head(consistent$Symbol, 50)
  } else {
    # Fallback: top genes by total significance
    all_cox |>
      filter(padj < COX_PADJ) |>
      count(Symbol, sort = TRUE) |>
      head(50) |>
      pull(Symbol)
  }

  if (length(top_genes_hm) >= 3) {
    hr_mat <- hr_wide |>
      filter(Symbol %in% top_genes_hm) |>
      column_to_rownames("Symbol") |>
      as.matrix()

    # Only keep cancers with results
    avail_cancers <- intersect(names(cox_per_cancer), colnames(hr_mat))
    hr_mat <- hr_mat[, avail_cancers, drop = FALSE]

    save_hr_heatmap(hr_mat, file.path(out_fig, "pancancer_hr_heatmap.png"))
    message("HR heatmap saved")
  }
}

# ============================================================================
# 7. Volcano-style plot: HR vs -log10(p) across all cancers
# ============================================================================

p_volcano <- all_cox |>
  filter(!is.na(HR), !is.na(pvalue)) |>
  mutate(
    log2HR    = log2(HR),
    log10p    = -log10(padj),
    Sig       = padj < COX_PADJ,
    Direction = ifelse(HR > 1, "High risk", "Protective")
  ) |>
  ggplot(aes(x = log2HR, y = log10p,
             color = interaction(Sig, Direction, sep = "_"))) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = -log10(COX_PADJ), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
  scale_color_manual(
    values = c(
      "FALSE_High risk"   = "gray70",
      "FALSE_Protective"  = "gray70",
      "TRUE_High risk"    = "#de2d26",
      "TRUE_Protective"   = "#2171b5"
    ),
    labels = c(
      "FALSE_High risk"   = "NS",
      "FALSE_Protective"  = "NS",
      "TRUE_High risk"    = "High risk (sig)",
      "TRUE_Protective"   = "Protective (sig)"
    ),
    name = NULL
  ) +
  facet_wrap(~ cancer, nrow = 1) +
  labs(
    title = "Pan-Cancer Cox Volcano Plot (Hub Genes)",
    x     = "Log2(Hazard Ratio)",
    y     = "-Log10(FDR)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    strip.text   = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(out_fig, "pancancer_cox_volcano.png"),
       plot = p_volcano,
       width  = 4 * length(cox_per_cancer),
       height = 5,
       dpi    = 300)

message("\nPan-cancer survival analysis complete.")
