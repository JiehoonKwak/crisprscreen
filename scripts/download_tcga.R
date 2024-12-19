library(TCGAbiolinks)

getGDCprojects()$project_id[grep("GBM", getGDCprojects()$project_id)]

proj_info <- getProjectSummary("TCGA-GBM")
proj_info[,c("file_count", "case_count", "data_category")]