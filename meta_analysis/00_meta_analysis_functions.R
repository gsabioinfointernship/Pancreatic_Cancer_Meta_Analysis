library(metafor)
library(metap)
library(tidyverse)
library(rio)

run_meta_analysis <- function(cancer, output_dir) {

  input_rds <- paste0("outputs/", cancer, "/meta_matrix/", cancer, "_meta_matrix_list.rds")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load and prepare list
  dataset_list <- readRDS(input_rds)

  diffexp <- lapply(dataset_list, function(df) {
    df |>
      filter(!is.na(log2FoldChange), !is.na(lfcSE), !is.na(pvalue),
             !is.na(Gene_Symbol), Gene_Symbol != "",
             lfcSE > 0, pvalue > 0, pvalue <= 1) |>
      select(Symbol = Gene_Symbol, Log2FC = log2FoldChange, SE = lfcSE, pvalue, padj) |>
      filter(!duplicated(Symbol))
  })

  all_genes <- unique(unlist(lapply(diffexp, function(df) df$Symbol)))
  message("  Total unique genes: ", length(all_genes))

  # Per-gene meta-analysis
  results <- lapply(seq_along(all_genes), function(i) {

    gene <- all_genes[i]
    if (i %% 1000 == 0) message("  Processing gene ", i, " / ", length(all_genes))

    # Extract gene data across studies
    gene_data <- lapply(diffexp, function(df) {
      row <- df[df$Symbol == gene, , drop = FALSE]
      if (nrow(row) == 0) return(NULL)
      row[1, ]
    })
    gene_data <- Filter(Negate(is.null), gene_data)
    n_studies  <- length(gene_data)
    if (n_studies < 2) return(NULL)

    yi    <- sapply(gene_data, function(x) x$Log2FC)
    sei   <- sapply(gene_data, function(x) x$SE)
    pvals <- sapply(gene_data, function(x) x$pvalue)

    # ── Random Effects Model
    # Try REML first; fall back to DL if Fisher scoring warns (boundary tau2=0)
    reml_warned <- FALSE
    rem <- withCallingHandlers(
      tryCatch(
        rma(yi = yi, sei = sei, method = "REML",
            control = list(maxiter = 1000, stepadj = 0.5)),
        error = function(e) NULL
      ),
      warning = function(w) {
        if (grepl("Fisher scoring|local maximum|tau\\^2", conditionMessage(w))) {
          reml_warned <<- TRUE
          invokeRestart("muffleWarning")
        }
      }
    )
    method_used <- "REML"
    # If REML warned and set tau2=0 with enough studies, try DL (closed-form, no iteration)
    if (reml_warned && !is.null(rem) && rem$tau2 == 0 && n_studies >= 4) {
      rem_dl <- tryCatch(
        suppressWarnings(rma(yi = yi, sei = sei, method = "DL")),
        error = function(e) NULL
      )
      if (!is.null(rem_dl) && rem_dl$tau2 > 0) {
        rem <- rem_dl
        method_used <- "DL_fallback"
      }
    }

    # ── Vote-Counting
    vote_up   <- sum(yi > 1  & pvals < 0.05, na.rm = TRUE)
    vote_down <- sum(yi < -1 & pvals < 0.05, na.rm = TRUE)

    # ── Fisher Combining ─────────────────────────────────────────────────────
    pvals_valid <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
    fisher_p <- if (length(pvals_valid) >= 2) {
      tryCatch(sumlog(pvals_valid)$p, error = function(e) NA_real_)
    } else NA_real_

    if (is.null(rem)) {
      data.frame(
        Symbol        = gene,
        n_studies     = n_studies,
        method_used   = NA_character_,
        pooled_log2FC = NA_real_,
        pooled_se     = NA_real_,
        rem_pvalue    = NA_real_,
        I2            = NA_real_,
        tau2          = NA_real_,
        Q             = NA_real_,
        Q_pval        = NA_real_,
        vote_up       = vote_up,
        vote_down     = vote_down,
        fisher_pvalue = fisher_p
      )
    } else {
      data.frame(
        Symbol        = gene,
        n_studies     = n_studies,
        method_used   = method_used,
        pooled_log2FC = as.numeric(rem$b),
        pooled_se     = as.numeric(rem$se),
        rem_pvalue    = as.numeric(rem$pval),
        I2            = rem$I2,
        tau2          = rem$tau2,
        Q             = rem$QE,
        Q_pval        = rem$QEp,
        vote_up       = vote_up,
        vote_down     = vote_down,
        fisher_pvalue = fisher_p
      )
    }
  })

  results_df <- bind_rows(Filter(Negate(is.null), results))

  # FDR correction
  results_df <- results_df |>
    mutate(
      rem_padj    = p.adjust(rem_pvalue,    method = "BH"),
      fisher_padj = p.adjust(fisher_pvalue, method = "BH")
    ) |>
    arrange(rem_pvalue)

  # Export
  export(results_df, paste0(output_dir, cancer, "_meta_analysis_results.csv"))
  saveRDS(results_df, paste0(output_dir, cancer, "_meta_analysis_results.rds"))
  message("  Saved: ", output_dir, cancer, "_meta_analysis_results.csv")

  invisible(results_df)
}
