# tabix -p bed fragments.tsv.gz 

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(AnnotationHub)
library(tidyverse)

# 1. Load Data
frags <- CreateFragmentObject(path = 'data/park2022_PO/ATAC/gbm/fragments.tsv.gz', cells = NULL)
peaks <- CallPeaks(frags, macs2.path = "macs2")

counts <- FeatureMatrix(fragments = frags, features = peaks, cells = NULL)
chrom_assay <- CreateChromatinAssay(counts = counts, fragments = fragment, peaks = peaks, sep = c("-", "-"), genome = "hg38")

obj <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")
obj

peaks_keep <- seqnames(granges(obj)) %in% standardChromosomes(granges(obj))
obj <- obj[as.vector(peaks_keep), ]

# 2. Annotation
ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens.v113")
ens <- ah[["AH119325"]]
annotations <- GetGrangesFromEnsDb(ensdb = ens)

saveRDS(obj, 'intermediate_output/park2022_PO_atac.rds')
