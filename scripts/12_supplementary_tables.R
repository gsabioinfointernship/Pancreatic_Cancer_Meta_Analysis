library(tidyverse)
library(rio)
library(writexl)

dir.create("outputs/SupplementaryTables", showWarnings = FALSE, recursive = TRUE)

supp_tables <- list()

geo_ids <- c("GSE130688", "GSE136569", "GSE171485", "GSE196009",
             "GSE211398", "GSE280271", "GSE293744")

dataset_summary <- tibble(
  GEO_Accession = geo_ids,
  Normal_Samples = c(15, 5, 6, 6, 12, 7, 8),
  Tumor_Samples = c(15, 5, 6, 13, 16, 6, 40),
  Total_Samples = c(30, 10, 12, 19, 28, 13, 48),
  Platform = c("Illumina NovaSeq 6000", "Illumina HiSeq 2500", "Illumina HiSeq 4000",
               "Illumina NovaSeq 6000", "Illumina NextSeq 500",
               "Illumina NextSeq 550", "Illumina NovaSeq 6000")
) %>%
  add_row(
    GEO_Accession = "Total",
    Normal_Samples = 59,
    Tumor_Samples = 101,
    Total_Samples = 160,
    Platform = "-"
  )

supp_tables[["Table_S1_Dataset_Characteristics"]] <- dataset_summary

if(file.exists("outputs/MetaVolcanoR/significant_genes.csv")) {

  sig_genes <- import("outputs/MetaVolcanoR/significant_genes.csv") %>%
    dplyr::select(Gene_ID, Gene_Symbol, randomSummary, randomP, randomP_fdr,
                  signcon, het_QEp) %>%
    dplyr::rename(
      Gene = Gene_ID,
      Symbol = Gene_Symbol,
      Log2FC_Meta = randomSummary,
      P_value = randomP,
      FDR = randomP_fdr,
      Sign_Consistency = signcon,
      Heterogeneity_P = het_QEp
    ) %>%
    arrange(P_value)

  supp_tables[["Table_S2_All_Significant_DEGs"]] <- sig_genes
}

if(file.exists("outputs/Enrichment/GO_BP_results.csv")) {

  go_results <- import("outputs/Enrichment/GO_BP_results.csv") %>%
    filter(p.adjust < 0.05) %>%
    dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
    dplyr::rename(
      GO_ID = ID,
      GO_Term = Description,
      Gene_Ratio = GeneRatio,
      Background_Ratio = BgRatio,
      P_value = pvalue,
      FDR = p.adjust,
      Gene_Count = Count
    ) %>%
    arrange(FDR)

  supp_tables[["Table_S3_GO_Enrichment"]] <- go_results
}

if(file.exists("outputs/Enrichment/KEGG_results.csv")) {

  kegg_results <- import("outputs/Enrichment/KEGG_results.csv") %>%
    filter(p.adjust < 0.05) %>%
    dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
    dplyr::rename(
      KEGG_ID = ID,
      Pathway = Description,
      Gene_Ratio = GeneRatio,
      Background_Ratio = BgRatio,
      P_value = pvalue,
      FDR = p.adjust,
      Gene_Count = Count
    ) %>%
    arrange(FDR)

  supp_tables[["Table_S4_KEGG_Pathways"]] <- kegg_results
}

if(file.exists("outputs/TCGA/validation_results.csv")) {

  tcga_validation <- import("outputs/TCGA/validation_results.csv") %>%
    dplyr::select(Gene_ID, Gene_Symbol, randomSummary, log2FoldChange,
                  padj, direction_concordance, validated) %>%
    dplyr::rename(
      Gene = Gene_ID,
      Symbol = Gene_Symbol,
      Meta_Log2FC = randomSummary,
      TCGA_Log2FC = log2FoldChange,
      TCGA_FDR = padj,
      Direction_Concordant = direction_concordance,
      Validated = validated
    ) %>%
    arrange(TCGA_FDR)

  supp_tables[["Table_S5_TCGA_Validation"]] <- tcga_validation
}

if(file.exists("outputs/Survival/Cox_univariate.csv")) {

  cox_univariate <- import("outputs/Survival/Cox_univariate.csv") %>%
    dplyr::select(Gene_ID, HR, HR_lower, HR_upper, p_value, FDR, concordance) %>%
    dplyr::rename(
      Gene = Gene_ID,
      Hazard_Ratio = HR,
      HR_95CI_Lower = HR_lower,
      HR_95CI_Upper = HR_upper,
      P_value = p_value,
      Concordance_Index = concordance
    ) %>%
    arrange(P_value)

  supp_tables[["Table_S6_Survival_Univariate"]] <- cox_univariate
}

if(file.exists("outputs/Survival/Cox_multivariate.csv")) {

  cox_multivariate <- import("outputs/Survival/Cox_multivariate.csv") %>%
    dplyr::select(Variable, HR, HR_lower, HR_upper, p_value) %>%
    dplyr::rename(
      Hazard_Ratio = HR,
      HR_95CI_Lower = HR_lower,
      HR_95CI_Upper = HR_upper,
      P_value = p_value
    )

  supp_tables[["Table_S7_Survival_Multivariate"]] <- cox_multivariate
}

if(file.exists("outputs/MachineLearning/model_performance.csv")) {

  ml_performance <- import("outputs/MachineLearning/model_performance.csv") %>%
    dplyr::select(Feature_Set, Model, Dataset, AUC, Accuracy,
                  Sensitivity, Specificity, PPV, NPV) %>%
    arrange(Feature_Set, Model, Dataset)

  supp_tables[["Table_S8_ML_Performance"]] <- ml_performance
}

if(file.exists("outputs/Network/hub_genes.csv")) {

  hub_genes <- import("outputs/Network/hub_genes.csv") %>%
    dplyr::select(Gene, Degree, Betweenness, Closeness, Eigenvector, Log2FC) %>%
    dplyr::rename(
      Degree_Centrality = Degree,
      Betweenness_Centrality = Betweenness,
      Closeness_Centrality = Closeness,
      Eigenvector_Centrality = Eigenvector,
      Meta_Log2FC = Log2FC
    ) %>%
    arrange(desc(Degree_Centrality))

  supp_tables[["Table_S9_Hub_Genes"]] <- hub_genes
}

if(file.exists("outputs/DrugTarget/druggable_genes.csv")) {

  druggable <- import("outputs/DrugTarget/druggable_genes.csv")

  if("n_drugs" %in% colnames(druggable)) {
    druggable_clean <- druggable %>%
      dplyr::select(geneName, n_drugs, drug_names, interaction_types, Log2FC, P_value, is_hub) %>%
      dplyr::rename(
        Gene = geneName,
        Number_of_Drugs = n_drugs,
        Drug_Names = drug_names,
        Interaction_Types = interaction_types,
        Meta_Log2FC = Log2FC,
        Hub_Gene = is_hub
      ) %>%
      arrange(desc(Number_of_Drugs))

  } else {
    druggable_clean <- druggable %>%
      dplyr::select(Gene_Symbol, Gene_Description, druggable_class, randomSummary, randomP, is_hub) %>%
      dplyr::rename(
        Gene = Gene_Symbol,
        Description = Gene_Description,
        Druggable_Class = druggable_class,
        Meta_Log2FC = randomSummary,
        P_value = randomP,
        Hub_Gene = is_hub
      ) %>%
      arrange(Druggable_Class, P_value)
  }

  supp_tables[["Table_S10_Druggable_Genes"]] <- druggable_clean
}

for(table_name in names(supp_tables)) {
  export(supp_tables[[table_name]],
         paste0("outputs/SupplementaryTables/", table_name, ".csv"))
}

write_xlsx(supp_tables, "outputs/SupplementaryTables/All_Supplementary_Tables.xlsx")

summary_table <- tibble(
  Table = names(supp_tables),
  Rows = sapply(supp_tables, nrow),
  Columns = sapply(supp_tables, ncol),
  Description = c(
    "Dataset characteristics for 7 GEO studies",
    "Complete list of 613 significant DEGs from meta-analysis",
    "GO Biological Process enrichment results (FDR < 0.05)",
    "KEGG pathway enrichment results (FDR < 0.05)",
    "TCGA-PAAD validation results for meta-analytical DEGs",
    "Univariate Cox regression survival analysis results",
    "Multivariate Cox regression results (genes + clinical variables)",
    "Machine learning model performance metrics",
    "Hub genes identified from PPI network analysis",
    "Druggable genes with drug interaction information"
  )[1:length(supp_tables)]
)

export(summary_table, "outputs/SupplementaryTables/Table_Summary.csv")

message("Supplementary tables created successfully!")
message("Total tables: ", length(supp_tables))
message("Output files:")
message("  - Individual CSV files in outputs/SupplementaryTables/")
message("  - Combined Excel file: All_Supplementary_Tables.xlsx")
