library(copykat)
library(Seurat)
library(here)
options(bitmapType = "cairo")

# 1. TNL
# obj <- readRDS(here("output/tnl.rds"))
# obj <- JoinLayers(obj)
# objs <- SplitObject(obj, split.by = "donor_id")

# ################
# for(i in seq_along(objs)) {
#     mat <- GetAssayData(objs[[i]], layer = "counts")

#     out <- copykat(rawmat = mat, sam.name = names(objs)[i], n.cores = 64)
#     pred <- data.frame(out$prediction)
#     predictions[[i]] <- pred
#     rm(mat, out)
#     gc()
# }
# combined_pred <- do.call(rbind, predictions)
# write.csv(combined_pred, here("output/copykat_tnl.csv"))


# 2. CPTAC
obj <- readRDS(here("output/cptac.rds"))
obj <- JoinLayers(obj)
objs <- SplitObject(obj, split.by = "donor_id")

predictions <- list()

for (i in seq_along(objs)) {
    mat <- GetAssayData(objs[[i]], layer = "counts")

    out <- copykat(rawmat = mat, sam.name = names(objs)[i], n.cores = 48)
    pred <- data.frame(out$prediction)
    predictions[[i]] <- pred
    rm(mat, out)
    gc()
}

combined_pred <- do.call(rbind, predictions)

write.csv(combined_pred, here("output/copykat_cptac.csv"))
