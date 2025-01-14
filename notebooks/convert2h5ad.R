library(Seurat)
library(sceasy)
library(reticulate)
library(here)
library(qs)
use_condaenv('/home/jiehoonk/micromamba/envs/sc_base', '/home/jiehoonk/.local/bin/micromamba') # TNL-A

# 1. TNL
obj <- qread(here("output/tnl2.qs"))
metadata <- obj@meta.data
obj[["mca"]] <- NULL
obj[["harmony"]] <- NULL
obj[["RNA"]] <- CreateAssayObject(counts = obj[["RNA"]]$counts)
obj@meta.data <- metadata
convertFormat(obj, from = "seurat", to = "anndata", outFile = here("output/tnl.h5ad"))

# 2. CPTAC
obj <- qread(here("output/cptac2.qs"))
metadata <- obj@meta.data
obj[["mca"]] <- NULL
obj[["harmony"]] <- NULL
obj[["RNA"]] <- CreateAssayObject(counts = obj[["RNA"]]$counts)
obj@meta.data <- metadata
convertFormat(obj, from = "seurat", to = "anndata", outFile = here("output/cptac.h5ad"))

# 3. WANG
obj <- qread(here("output/wang2024_2.qs"))
metadata <- obj@meta.data
obj[["mca"]] <- NULL
obj[["harmony"]] <- NULL
obj[["RNA"]] <- CreateAssayObject(counts = obj[["RNA"]]$counts)
obj@meta.data <- metadata
convertFormat(obj, from = "seurat", to = "anndata", outFile = here("output/wang2024.h5ad"))