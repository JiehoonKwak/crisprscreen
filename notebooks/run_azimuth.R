library(Seurat)
library(Azimuth)
library(here)
library(qs)
source(here('data/gbmap/modified_azimuth.R'))
ref <- readRDS(here("data/gbmap/azimuth_core_GBmap.rds"))

# TNL - DONE
obj <- readRDS(here("output/tnl2.rds"))
obj <- RunAzimuth(obj, reference = ref, annotation.levels = c('annotation_level_3', 'annotation_level_4'))

qsave(obj, here("output/tnl2.qs"))

# CPTAC - DONE
obj <- readRDS(here("output/cptac2.rds"))
obj <- RunAzimuth(obj, reference = ref, annotation.levels = c('annotation_level_3', 'annotation_level_4'))

qsave(obj, here("output/cptac2.qs"))


# WANG2024 - DONE
obj <- readRDS(here("output/wang2024_2.rds"))
obj <- RunAzimuth(obj, reference = ref, annotation.levels = c('annotation_level_3', 'annotation_level_4'))

qsave(obj, here("output/wang2024_2.qs"))

