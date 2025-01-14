library(Seurat)
library(Azimuth)
ref <- readRDS(here("data/gbmap/azimuth_core_GBmap.rds"))
obj <- RunAzimuth(obj, reference = ref, annotation.levels = c('annotation_level_3', 'annotation_level_4'))