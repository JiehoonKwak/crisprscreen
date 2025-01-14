library(Seurat)
library(tidyverse)
library(here)
library(SingleR)
library(celldex)
library(qs)

# TNL - DONE
obj <- readRDS(here('output/tnl2.rds'))
ref <- HumanPrimaryCellAtlasData()
counts <- GetAssayData(obj, layer = "counts")
pred <- SingleR(test = counts, ref = ref, labels = ref$label.main)

obj <- AddMetaData(obj, pred$pruned.labels, col.name = "SingleR")

saveRDS(obj, here("output/tnl2.rds"))


# CPTAC - DONE
obj <- readRDS(here('output/cptac2.rds'))
ref <- HumanPrimaryCellAtlasData()
counts <- GetAssayData(obj, layer = "counts")
pred <- SingleR(test = counts, ref = ref, labels = ref$label.main)

obj <- AddMetaData(obj, pred$pruned.labels, col.name = "SingleR")

saveRDS(obj, here("output/cptac2.rds"))


# WANG2024
obj <- qread(here('output/wang2024_2.qs'))
ref <- HumanPrimaryCellAtlasData()
counts <- GetAssayData(obj, layer = "counts")
pred <- SingleR(test = counts, ref = ref, labels = ref$label.main)

obj <- AddMetaData(obj, pred$pruned.labels, col.name = "SingleR")

obj[[]]

qsave(obj, here("output/wang2024_2.qs")) #send to b
