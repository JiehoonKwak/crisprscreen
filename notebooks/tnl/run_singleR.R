library(Seurat)
library(tidyverse)
library(here)
library(SingleR)
library(celldex)

# TNL - DONE
obj <- readRDS(here('output/tnl2.rds'))
ref <- HumanPrimaryCellAtlasData()
counts <- GetAssayData(obj, layer = "counts")
pred <- SingleR(test = counts, ref = ref, labels = ref$label.main)

obj <- AddMetaData(obj, pred$pruned.labels, col.name = "SingleR")

saveRDS(obj, here("output/tnl2.rds"))


# CPTAC
obj <- readRDS(here('output/cptac2.rds'))
ref <- HumanPrimaryCellAtlasData()
counts <- GetAssayData(obj, layer = "counts")
pred <- SingleR(test = counts, ref = ref, labels = ref$label.main)

obj <- AddMetaData(obj, pred$pruned.labels, col.name = "SingleR")

saveRDS(obj, here("output/cptac2.rds"))

obj[[]]

# WANG2024
obj <- readRDS(here('output/wang2024_2.rds'))
ref <- HumanPrimaryCellAtlasData()
counts <- GetAssayData(obj, layer = "counts")
pred <- SingleR(test = counts, ref = ref, labels = ref$label.main)

obj <- AddMetaData(obj, pred$pruned.labels, col.name = "SingleR")

saveRDS(obj, here("output/wang2024_2.rds"))