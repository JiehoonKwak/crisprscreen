library(copykat)
library(Seurat)
library(here)

obj <- readRDS(here("output/tnl.rds"))
obj <- JoinLayers(obj)
objs <- SplitObject(obj, split.by = "donor_id")

################
for(i in seq_along(objs)) {
    mat <- GetAssayData(objs[[i]], layer = "counts")
    
    out <- copykat(rawmat = mat, sam.name = names(objs)[i], n.cores = 64)
    pred <- data.frame(out$prediction)
    predictions[[i]] <- pred
    rm(mat, out)
    gc()
}
combined_pred <- do.call(rbind, predictions)
write.csv(combined_pred, here("output/copykat_tnl.csv"))








##################33
mats <- lapply(objs, function(obj) {
    mat <- GetAssayData(obj, layer = "counts")
    return(mat)
})

outs <- lapply(mats, function(mat) {
    copykat(rawmat = mat, sam.name = names(mat), n.cores = 64)
})
predictions <- lapply(outs, function(out) data.frame(out$prediction))

objs <- mapply(function(obj, out) {
    pred.test <- data.frame(out$prediction)
    AddMetaData(obj, pred.test, col.name = "copykat")
}, objs, outs, SIMPLIFY = FALSE)

names(objs)

tnl <- merge(objs[[1]], y = unlist(objs[2:length(objs)]), add.cell.ids = names(objs), project = "tnl", merge.data = T)