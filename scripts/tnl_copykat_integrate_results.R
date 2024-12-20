library(tidyverse)
output_dir <- "/home/jiehoonk/project/crisprscreen/intermediate_output/copykat_results"
files <- list.files(output_dir, full.names = TRUE)

obj <- readRDS(files[1])
cnv_data <- obj$CNAmat %>% as.data.frame()
cnv_scores <- colMeans(cnv_data, na.rm = TRUE) 
cnv_scores  |>  select(str_starts(.$names, "sample_")) |> as.data.frame() 
