# GBM TF Analysis Pipeline

# Load required packages
library(tidyverse)
library(Seurat)
library(biomaRt)
library(DESeq2)
library(fgsea)
library(TFBSTools)
library(JASPAR2020)

# 1. Data Loading and Processing
load_geo_data <- function(geo_accession) {
  # Load GEO data using GEOquery
  gset <- getGEO(geo_accession, GSEMatrix = TRUE)
  # Extract expression matrix
  expr_mat <- exprs(gset[[1]])
  # Get metadata
  metadata <- pData(gset[[1]])
  return(list(expr = expr_mat, meta = metadata))
}

# 2. TF List Generation
get_tf_list <- function() {
  # Get human TFs from AnimalTFDB
  human_tfs <- read.table("http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF.txt",
                         header = TRUE, sep = "\t")
  return(human_tfs$Symbol)
}

# 3. Differential Expression Analysis
run_de_analysis <- function(counts, metadata, tf_list) {
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = metadata,
                               design = ~ condition)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds)
  
  # Filter for TFs
  tf_results <- res[rownames(res) %in% tf_list,]
  
  return(tf_results)
}

# 4. Integration with Normal Brain Data
compare_with_normal <- function(tumor_expr, normal_expr, tf_list) {
  # Normalize data
  tumor_norm <- normalize_data(tumor_expr)
  normal_norm <- normalize_data(normal_expr)
  
  # Compare expression
  tf_diff <- calculate_tf_differences(tumor_norm, normal_norm, tf_list)
  
  return(tf_diff)
}

# 5. Survival Analysis
survival_analysis <- function(expr_data, clinical_data, tf_list) {
  library(survival)
  library(survminer)
  
  # Perform survival analysis for each TF
  tf_survival <- lapply(tf_list, function(tf) {
    fit <- coxph(Surv(time, status) ~ expr_data[tf,], data = clinical_data)
    return(summary(fit))
  })
  
  return(tf_survival)
}

# 6. TF Ranking and Prioritization
rank_tfs <- function(de_results, survival_results, conservation_scores) {
  # Combine multiple metrics
  tf_scores <- data.frame(
    tf = rownames(de_results),
    de_score = -log10(de_results$padj) * sign(de_results$log2FoldChange),
    survival_score = unlist(lapply(survival_results, function(x) x$conf.int[1])),
    conservation = conservation_scores
  )
  
  # Calculate composite score
  tf_scores$composite_score <- scale(tf_scores$de_score) +
    scale(tf_scores$survival_score) +
    scale(tf_scores$conservation)
  
  # Rank TFs
  ranked_tfs <- tf_scores[order(-tf_scores$composite_score),]
  
  return(ranked_tfs)
}

# 7. Export Results
export_results <- function(ranked_tfs, output_dir) {
  # Export ranked TF list
  write_csv(ranked_tfs, file.path(output_dir, "ranked_tfs.csv"))
  
  # Generate summary plots
  plot_tf_rankings(ranked_tfs)
  
  # Export gRNA suggestions
  export_grna_suggestions(ranked_tfs$tf[1:100])
}

# Main workflow
main <- function() {
  # Load and process data
  gbm_data <- load_geo_data("GSE67835")
  normal_data <- load_geo_data("GSE99333")
  
  # Get TF list
  tf_list <- get_tf_list()
  
  # Run differential expression
  de_results <- run_de_analysis(gbm_data$expr, gbm_data$meta, tf_list)
  
  # Compare with normal tissue
  normal_comparison <- compare_with_normal(gbm_data$expr, normal_data$expr, tf_list)
  
  # Run survival analysis
  surv_results <- survival_analysis(gbm_data$expr, gbm_data$meta, tf_list)
  
  # Get conservation scores
  conservation_scores <- get_conservation_scores(tf_list)
  
  # Rank TFs
  ranked_tfs <- rank_tfs(de_results, surv_results, conservation_scores)
  
  # Export results
  export_results(ranked_tfs, "output")
}