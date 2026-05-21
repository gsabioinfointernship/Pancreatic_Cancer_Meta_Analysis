library(tidyverse)
library(rio)
library(httr)
library(jsonlite)

select <- dplyr::select
filter <- dplyr::filter

CANCERS         <- c("colorectal", "esophagus", "kidney", "liver", "pancreatic")
PADJ_CUTOFF     <- 0.05
LFC_CUTOFF      <- 1
MIN_CANCERS_PC  <- 3

# ── Druggable gene family keywords
DRUGGABLE_FAMILIES <- c(
  "kinase", "receptor", "channel", "transporter", "protease",
  "phosphatase", "transferase", "hydrolase", "oxidoreductase",
  "GPCR", "nuclear receptor", "integrin", "ATPase", "helicase",
  "methyltransferase", "dehydrogenase", "ligase", "polymerase"
)

# ── Load meta-analysis results
load_meta <- function(cancer) {
  path <- paste0("outputs/", cancer, "/meta_analysis/", cancer, "_meta_analysis_results.rds")
  readRDS(path)
}

# ── Load hub genes (returns NULL if network not yet run)
load_hub_genes <- function(cancer) {
  path <- paste0("outputs/", cancer, "/network/", cancer, "_hub_genes.csv")
  if (!file.exists(path)) return(NULL)
  import(path)
}

# ── Query DGIdb v5 GraphQL API
# Batches to keep query size manageable
query_dgidb <- function(gene_symbols, batch_size = 50) {
  if (length(gene_symbols) == 0) return(NULL)

  gene_symbols <- unique(gene_symbols[!is.na(gene_symbols) & gene_symbols != ""])
  batches      <- split(gene_symbols, ceiling(seq_along(gene_symbols) / batch_size))
  url          <- "https://dgidb.org/api/graphql"

  results <- lapply(seq_along(batches), function(i) {
    batch <- batches[[i]]
    message("    DGIdb batch ", i, "/", length(batches),
            " (", length(batch), " genes)...")

    # Build GraphQL query
    gene_list <- paste0('"', batch, '"', collapse = ", ")
    query <- paste0(
      '{ genes(names: [', gene_list, ']) {',
      '  nodes {',
      '    name',
      '    interactions {',
      '      drug { name approved }',
      '      interactionScore',
      '      interactionTypes { type directionality }',
      '      sources { sourceDbName }',
      '      publications { pmid }',
      '    }',
      '  }',
      '} }'
    )

    resp <- tryCatch(
      POST(
        url,
        body    = toJSON(list(query = query), auto_unbox = TRUE),
        add_headers("Content-Type" = "application/json"),
        timeout(90)
      ),
      error = function(e) { message("    HTTP error: ", e$message); NULL }
    )

    if (is.null(resp) || http_error(resp)) {
      message("    DGIdb request failed (status ",
              if (!is.null(resp)) status_code(resp) else "NA", ")")
      return(NULL)
    }

    parsed <- tryCatch(
      fromJSON(content(resp, "text", encoding = "UTF-8")),
      error = function(e) { message("    JSON parse error: ", e$message); NULL }
    )
    if (is.null(parsed)) return(NULL)

    nodes <- tryCatch(parsed$data$genes$nodes, error = function(e) NULL)
    if (is.null(nodes) || length(nodes) == 0) return(NULL)

    # Unnest per gene
    rows <- lapply(seq_len(nrow(nodes)), function(j) {
      symbol <- nodes$name[j]
      iacts  <- nodes$interactions[[j]]
      if (is.null(iacts) || length(iacts) == 0 || !is.data.frame(iacts) || nrow(iacts) == 0) return(NULL)

      tibble(
        Symbol          = symbol,
        Drug            = iacts$drug$name,
        Approved        = iacts$drug$approved,
        InteractionType = sapply(iacts$interactionTypes,
                                 function(x) paste(unlist(x[["type"]]), collapse = "; ")),
        Sources         = sapply(iacts$sources,
                                 function(x) paste(unlist(x[["sourceDbName"]]), collapse = "; ")),
        Score           = iacts$interactionScore,
        PMIDs           = sapply(iacts$publications,
                                 function(x) paste(unlist(x[["pmid"]]), collapse = "; "))
      )
    })
    bind_rows(rows)
  })

  out <- bind_rows(results)
  if (nrow(out) == 0) return(NULL)
  out
}

# ── Classify drug approval status from source/drug name patterns
classify_approval <- function(drug_name, sources = "") {
  drug_lower   <- tolower(drug_name)
  source_lower <- tolower(sources)
  if (str_detect(source_lower, "fda|drugbank approved|chembl approved")) return("Approved")
  if (str_detect(drug_lower, "^[A-Z]{3,}-\\d+$|phase|investig")) return("Investigational")
  "Unknown"
}

# ── Classify druggable gene family by symbol
classify_gene_family <- function(symbol) {
  s <- toupper(symbol)
  case_when(
    str_detect(s, "^CDK|^MAPK|^EGFR|^ERBB|^MET$|^ALK$|^RET$|^FGFR|^PDGFR|^VEGFR|KIT$|^ABL|^JAK|^SRC$|^BTK$|^PIK3|^AKT|^MTOR$|^BRAF$|^RAF|^ROS1$|^NTRK|KLK|KINASE") ~ "Kinase",
    str_detect(s, "^MMP|PROT[EA]|^CASP|TRYPSIN|ELASTASE|PRSS|CTSL|CTSB|FURIN") ~ "Protease",
    str_detect(s, "^VEGFA?$|^IGF|^TGFB|^TNF$|^IL[0-9]|^CSF|^EGF$|^FGF|^WNT|^HGF$|^PDGF") ~ "Growth Factor/Cytokine",
    str_detect(s, "RECEPTOR|^EGFR$|^ERBB|^ESR|^AR$|^RXRA?|^PPARG|^NR[0-9]|^PTCH|^FZD") ~ "Receptor",
    str_detect(s, "^ABCB|^ABCC|^ABCG|^SLC[0-9]|TRANSPORTER|^MDR|^CFTR") ~ "Transporter",
    str_detect(s, "^KCNQ|^SCN|^CACNA|^HCN|CHANNEL") ~ "Ion Channel",
    str_detect(s, "^DNMT|^EZH|^HDAC|^KDM|^KMT|^PRMT|^SIRT|EPIGENETIC|METHYLTRANSFER|HISTONE") ~ "Epigenetic",
    str_detect(s, "^PARP|^TOP[12]|^TERT$|^POLR|POLYMERASE|HELICASE") ~ "DNA Repair/Replication",
    str_detect(s, "^MDM2?$|^BCL2$|^BCL[X2]|^MCL1$|^XIAP$|APOPTOSIS") ~ "Apoptosis Regulator",
    str_detect(s, "^KRAS$|^HRAS$|^NRAS$|^RHOA?$|^RAC1$|^CDC42$|GTPASE") ~ "GTPase/RAS",
    str_detect(s, "^HIF1|^VEGF|ANGIO") ~ "Angiogenesis",
    TRUE ~ "Other"
  )
}

# ── Build prioritised drug target table
build_priority_table <- function(iact_df, meta_df, hub_df = NULL) {
  if (is.null(iact_df) || nrow(iact_df) == 0) return(NULL)

  # Per-gene drug summary
  by_gene <- iact_df |>
    group_by(Symbol) |>
    summarise(
      n_drugs          = n_distinct(Drug),
      n_interactions   = n(),
      drug_names       = paste(unique(Drug), collapse = "; "),
      interaction_types = paste(unique(InteractionType), collapse = "; "),
      sources          = paste(unique(Sources), collapse = "; "),
      .groups = "drop"
    )

  # Join meta data
  meta_slim <- meta_df |>
    select(Symbol, pooled_log2FC, rem_pvalue, rem_padj) |>
    filter(!is.na(Symbol))

  priority <- by_gene |>
    left_join(meta_slim, by = "Symbol") |>
    mutate(
      Direction    = ifelse(pooled_log2FC > 0, "Up", "Down"),
      Gene_Family  = classify_gene_family(Symbol),
      is_hub       = if (!is.null(hub_df)) Symbol %in% hub_df$Symbol else FALSE,
      Priority     = case_when(
        is_hub & n_drugs >= 3 ~ "High",
        is_hub | n_drugs >= 3 ~ "Medium",
        TRUE                  ~ "Low"
      )
    ) |>
    arrange(desc(is_hub), desc(n_drugs))

  priority
}

# ── Save drug count lollipop plot (top N genes)
save_drug_lollipop <- function(priority_df, title, path, n = 25) {
  if (is.null(priority_df) || nrow(priority_df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- head(priority_df, n)
  p <- ggplot(df, aes(x = reorder(Symbol, n_drugs), y = n_drugs, color = Direction)) +
    geom_segment(aes(xend = Symbol, y = 0, yend = n_drugs), linewidth = 0.8) +
    geom_point(aes(size = ifelse(is_hub, 4, 2.5)), show.legend = FALSE) +
    geom_point(aes(size = ifelse(is_hub, 4, 2.5), color = Direction)) +
    scale_color_manual(values = c(Up = "#de2d26", Down = "#2171b5"), name = "Direction") +
    scale_size_identity() +
    coord_flip() +
    labs(title = title, x = NULL, y = "Number of Drugs") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(path, plot = p, width = 9, height = 7, dpi = 300)
  invisible(p)
}

# ── Save gene family barplot
save_family_barplot <- function(priority_df, title, path) {
  if (is.null(priority_df) || nrow(priority_df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- priority_df |>
    count(Gene_Family, name = "n_genes") |>
    filter(Gene_Family != "Other") |>
    arrange(desc(n_genes))
  if (nrow(df) == 0) return(invisible(NULL))
  p <- ggplot(df, aes(x = reorder(Gene_Family, n_genes), y = n_genes, fill = Gene_Family)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_text(aes(label = n_genes), hjust = -0.2, size = 3.5) +
    coord_flip() +
    labs(title = title, x = NULL, y = "Number of Druggable Genes") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(path, plot = p, width = 8, height = 6, dpi = 300)
  invisible(p)
}

# ── Save hub vs non-hub drug count boxplot
save_hub_drug_boxplot <- function(priority_df, title, path) {
  if (is.null(priority_df) || nrow(priority_df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  p <- ggplot(priority_df, aes(x = is_hub, y = n_drugs, fill = is_hub)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
    scale_fill_manual(values = c("FALSE" = "#bdbdbd", "TRUE" = "#fd8d3c"),
                      labels = c("Non-Hub", "Hub"), name = NULL) +
    scale_x_discrete(labels = c("FALSE" = "Non-Hub", "TRUE" = "Hub")) +
    labs(title = title, x = NULL, y = "Number of Drugs") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none")
  ggsave(path, plot = p, width = 6, height = 5, dpi = 300)
  invisible(p)
}

# ── Safe CSV export: converts ALL list columns to character, bypasses fwrite
safe_export <- function(df, path) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- sapply(df[[nm]], function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        tryCatch(
          paste(unique(rapply(x, as.character, how = "unlist")), collapse = "; "),
          error = function(e) as.character(x[[1]])[1]
        )
      }, USE.NAMES = FALSE)
    }
  }
  write.csv(df, path, row.names = FALSE)
}

# ── Export outputs
export_drug <- function(iact_df, priority_df, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  safe_export(iact_df,     file.path(out_dir, paste0(prefix, "_drug_interactions.csv")))
  safe_export(priority_df, file.path(out_dir, paste0(prefix, "_drug_targets.csv")))
}
