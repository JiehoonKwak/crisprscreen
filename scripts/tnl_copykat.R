#!/usr/bin/env Rscript

# Load necessary libraries
library(schard)
library(here)
library(Seurat)
library(copykat)

# Define the input file
input_file <- here("output/annotated.h5ad")
print(input_file)

# Convert the input file to a Seurat object
obj <- schard::h5ad2seurat(input_file)

# Split the Seurat object based on sample_id
sample_id <- obj@meta.data$sample_id
seurat_list <- SplitObject(obj, split.by = "sample_id")

# Output directory
output_dir <- "/home/jiehoonk/project/crisprscreen/intermediate_output/copykat_results"
dir.create(output_dir, showWarnings = FALSE)

# Iterate over each subset and run the analysis
for (sample_name in names(seurat_list)) {
  message("Processing sample: ", sample_name)
  
  # Get the subset
  subset_obj <- seurat_list[[sample_name]]
  
  # Extract the raw data matrix
  raw <- as.matrix(LayerData(subset_obj, layer = "data"))
  
  # Run copykat
  copykat_result <- copykat(rawmat = raw, genome = "hg20", n.cores = 48)
  
  # Save results
  result_file <- file.path(output_dir, paste0("copykat_result_", sample_name, ".rds"))
  saveRDS(copykat_result, file = result_file)
  
  message("Results saved to: ", result_file)
}

message("All samples processed.")