library(Seurat)
library(tidyverse)

list.files('data/park2022_PO/RNA')

rna <- Read10X('data/park2022_PO/RNA/primaryGBM')
str(rna)
